#include "ILP_cplex.h"

std::vector<int> solveWithCplexST(const CplexProblemData& data, 
                                double mipGap, double timeLimit,
                                const std::vector<int>& initialRoute,                                
                                double& obj_bound, double& result_gap) {
    std::vector<int> route;
    try {
        IloEnv env;
        IloModel model(env);
        int N = data.N;

        IloArray<IloBoolVarArray> x(env, N);
        for (int i = 0; i < N; ++i)
            x[i] = IloBoolVarArray(env, N);

        IloNumVarArray u(env, N, 0, N, ILOINT);

        IloExpr obj(env);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j)
                    obj += data.Prize[j] * x[i][j];
            }
        }
        model.add(IloMaximize(env, obj));
        obj.end();

        // At_Most_One_Incoming and At_Most_One_Outgoing
        for (int i = 0; i < N; ++i) {
            IloExpr in(env), out(env);
            for (int j = 0; j < N; ++j) {
                if (j != i) {
                    in += x[j][i];
                    out += x[i][j];
                }
            }
            model.add(in <= 1);
            model.add(out <= 1);
            in.end(); out.end();
        }

        // Flow_Balance
        IloExpr startOut(env);
        for (int j = 0; j < N; ++j) {
            if (j != data.startIndex)
                startOut += x[data.startIndex][j];
        }
        model.add(startOut == 1);
        startOut.end();

        for (int i = 0; i < N; ++i) {
            if (i != data.startIndex)
                model.add(x[i][data.startIndex] == 0);
        }

        IloExpr endIn(env);
        for (int i = 0; i < N; ++i) {
            if (i != data.endIndex)
                endIn += x[i][data.endIndex];
        }
        model.add(endIn == 1);
        endIn.end();

        for (int j = 0; j < N; ++j) {
            if (j != data.endIndex)
                model.add(x[data.endIndex][j] == 0);
        }

        for (int i = 0; i < N; ++i) {
            if (i == data.startIndex || i == data.endIndex) continue;
            IloExpr flow(env);
            for (int j = 0; j < N; ++j) {
                if (j != i) flow += x[j][i] - x[i][j];
            }
            model.add(flow == 0);
            flow.end();
        }

        // Budget Constraint
        IloExpr totalCost(env);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j)
                    totalCost += data.Cost[i][j] * x[i][j];
            }
        }
        model.add(totalCost <= data.Budget);
        totalCost.end();

        // Subtour_Elimination (MTZ)
        for (int i = 0; i < N; ++i) {
            if (i == data.startIndex) continue;
            for (int j = 0; j < N; ++j) {
                if (j == data.startIndex || i == j) continue;
                model.add(u[i] - u[j] + N * x[i][j] <= N - 1);
            }
        }

        IloCplex cplex(model);
        if (mipGap >= 0)
            cplex.setParam(IloCplex::EpGap, mipGap);
        if (timeLimit >= 0)
            cplex.setParam(IloCplex::TiLim, timeLimit);

        // Add initial solution if provided
        if (!initialRoute.empty()) {
            IloNumVarArray startVars(env);
            IloNumArray startVals(env);

            for (size_t idx = 0; idx + 1 < initialRoute.size(); ++idx) {
                int i = initialRoute[idx];
                int j = initialRoute[idx + 1];
                startVars.add(x[i][j]);
                startVals.add(1.0);
            }

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    bool inRoute = false;
                    for (size_t idx = 0; idx + 1 < initialRoute.size(); ++idx) {
                        if (initialRoute[idx] == i && initialRoute[idx + 1] == j) {
                            inRoute = true;
                            break;
                        }
                    }
                    if (!inRoute) {
                        startVars.add(x[i][j]);
                        startVals.add(0.0);
                    }
                }
            }

            for (size_t k = 1; k + 1 < initialRoute.size(); ++k) {
                int i = initialRoute[k];
                startVars.add(u[i]);
                startVals.add(k + 1);
                startVals.add(k + 2);
            }

            cplex.addMIPStart(startVars, startVals);
            startVars.end(); startVals.end();
        }

        if (!cplex.solve()) {
            env.error() << "Failed to optimize." << std::endl;
            throw(-1);
        }

        if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible) {
            obj_bound = cplex.getBestObjValue();
            result_gap = cplex.getMIPRelativeGap();
            std::cout << "Best bound: " << obj_bound << " (Status: " << cplex.getStatus() << ")\n";
            std::cout << "Optimality gap: " << result_gap << std::endl;
        } else {
            std::cout << "No optimal or feasible solution found.\n";
        }

        env.out() << "Objective value: " << cplex.getObjValue() << std::endl;

        int current = data.startIndex;
        route.push_back(current);
        while (current != data.endIndex) {
            bool found = false;
            for (int j = 0; j < N; ++j) {
                if (j != current && cplex.getValue(x[current][j]) > 0.5) {
                    route.push_back(j);
                    current = j;
                    found = true;
                    break;
                }
            }
            if (!found) break;
        }

        env.end();
    } catch (IloException& e) {
        std::cerr << "Cplex exception: " << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error occurred in CPLEX solver." << std::endl;
    }

    return route;
}


std::vector<int> cplexSolveST(const std::vector<std::vector<double>>& Cost,
                               const std::vector<double>& Prize,
                               int start, int end, double Budget,
                               double& obj_value, double& result_gap,
                               double mipGap, double timeLimit,
                               const std::vector<int>& initialRoute) {
    int N = Prize.size();
    CplexProblemData data(N, start, end, Budget, Cost, Prize);
    return solveWithCplexST(data, mipGap, timeLimit, initialRoute, obj_value, result_gap);
}
