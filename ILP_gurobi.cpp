#include "ILP_gurobi.h"
#include <string>
#include <iostream>
#include <cmath>
#include <limits>

std::vector<int> solveWithGurobi(const ProblemData& data, double mipGap, double timeLimit,
                                 const std::vector<int>& initialRoute) {
    std::vector<int> route;
    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_tsp.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TSP_PrizeCollecting");

        // Set solver parameters
        // model.set(GRB_DoubleParam_MIPGap, 0.001); 
        // model.set(GRB_DoubleParam_TimeLimit, 500.0); 
        std::cout << "mipGap: " << mipGap << "timeLimit" << timeLimit << std::endl;
        if (mipGap != -1)
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit != -1)
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);

        int N = data.N;
        int start = data.startIndex;
        double Budget = data.Budget;
        const auto& Cost = data.Cost;
        const auto& Prize = data.Prize;

        // x[i][j]
        std::vector<std::vector<GRBVar>> x(N, std::vector<GRBVar>(N));
        // is_end[i]
        std::vector<GRBVar> is_end(N);
        // u[i]
        std::vector<GRBVar> u(N);

        // 1) x[i][j] with objective = Prize[j] if i != j, else 0
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                double objCoeff = (i == j) ? 0.0 : Prize[j];
                x[i][j] = model.addVar(0.0, 1.0, objCoeff, GRB_BINARY,
                                       "x_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // 2) is_end[i]
        for(int i = 0; i < N; i++) {
            is_end[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     "is_end_" + std::to_string(i));
        }

        // 3) u[i]
        for(int i = 0; i < N; i++) {
            if(i == start) continue;
            u[i] = model.addVar(2.0, (double)N, 0.0, GRB_INTEGER,
                                "u_" + std::to_string(i));
        }
        if (!initialRoute.empty()) {
            if (initialRoute[0] != start) {
                std::cerr << "Initial route must start at " << start << "\n";
                return route;
            }
            double totalCost = 0.0;
            for (size_t k = 0; k < initialRoute.size() - 1; k++) {
                totalCost += Cost[initialRoute[k]][initialRoute[k + 1]];
            }
            if (totalCost > Budget) {
                std::cerr << "Initial route exceeds budget: " << totalCost << " > " << Budget << "\n";
                return route;
            }
            // Reset all x[i][j]
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    x[i][j].set(GRB_DoubleAttr_Start, 0.0);
                }
                is_end[i].set(GRB_DoubleAttr_Start, 0.0);
            }
            // Set initial route
            for (size_t k = 0; k < initialRoute.size() - 1; k++) {
                int i = initialRoute[k];
                int j = initialRoute[k + 1];
                x[i][j].set(GRB_DoubleAttr_Start, 1.0);
            }
            int endNode = initialRoute.back();
            is_end[endNode].set(GRB_DoubleAttr_Start, 1.0);
        }

        // 5) Objective MAXIMIZE
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        // Constraints
        // A) One_End_Node
        {
            GRBLinExpr lhs = 0.0;
            for(int i = 0; i < N; i++) {
                lhs += is_end[i];
            }
            model.addConstr(lhs == 1, "One_End_Node");
        }
        // B) At_Most_One_Incoming
        for(int i = 0; i < N; i++) {
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++) {
                if(j != i) lhs += x[j][i];
            }
            model.addConstr(lhs <= 1.0, "AtMostOneIncoming_" + std::to_string(i));
        }
        // C) At_Most_One_Outgoing
        for(int i = 0; i < N; i++) {
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++) {
                if(j != i) lhs += x[i][j];
            }
            model.addConstr(lhs <= 1.0, "AtMostOneOutgoing_" + std::to_string(i));
        }
        // D) Flow_Balance
        for(int i = 0; i < N; i++) {
            if(i == start) continue;
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++) {
                if(j != i) lhs += x[i][j];
            }
            for(int j = 0; j < N; j++) {
                if(j != i) lhs -= x[j][i];
            }
            lhs += is_end[i];
            model.addConstr(lhs == 0.0, "FlowBalance_" + std::to_string(i));
        }
        // E) Start_Flow
        {
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++) {
                if(j != start) lhs += x[start][j];
            }
            for(int j = 0; j < N; j++) {
                if(j != start) lhs -= x[j][start];
            }
            model.addConstr(lhs == 1.0, "StartFlow");
        }
        // F) Budget_Constraint
        {
            GRBLinExpr lhs = 0.0;
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++) {
                    if(i != j) {
                        lhs += Cost[i][j] * x[i][j];
                    }
                }
            }
            model.addConstr(lhs <= Budget, "BudgetConstraint");
        }
        // G) Subtour_Elimination (MTZ)
        for(int i = 0; i < N; i++) {
            if(i == start) continue;
            for(int j = 0; j < N; j++) {
                if(j == i || j == start) continue;
                GRBLinExpr lhs = 0.0;
                lhs += u[i];
                lhs -= u[j];
                lhs += (double)N * x[i][j];
                model.addConstr(lhs <= (double)N - 1.0,
                                "MTZ_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // Optimize with the initial solution
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if(status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            double objVal = model.get(GRB_DoubleAttr_ObjVal);
            std::cout << "Objective value: " << objVal << " (Status: " << status << ")\n";

            // Build the route
            route.push_back(start);
            int currentNode = start;
            while(true) {
                bool foundNext = false;
                for(int nextNode = 0; nextNode < N; nextNode++) {
                    if(nextNode == currentNode) continue;
                    double val = x[currentNode][nextNode].get(GRB_DoubleAttr_X);
                    if(val > 0.5) {
                        route.push_back(nextNode);
                        currentNode = nextNode;
                        foundNext = true;
                        break;
                    }
                }
                if(!foundNext) {
                    break;
                }
            }
        } else {
            std::cout << "No optimal solution found. Status = " << status << "\n";
        }

    } catch(GRBException &e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
        return route;
    } catch(...) {
        std::cerr << "Unknown error in solveWithGurobi.\n";
        return route;
    }

    return route;
}

std::vector<int> gurobiSolve(const std::vector<std::vector<double>>& Cost, 
                             const std::vector<double>& Prize,
                             int start, double Budget, double mipGap, double timeLimit,
                             const std::vector<int>& initialRoute) {
    int N = Prize.size();
    ProblemData data(N, start, Budget, Cost, Prize);
    std::vector<int> path = solveWithGurobi(data, mipGap, timeLimit, initialRoute);
    return path;
}