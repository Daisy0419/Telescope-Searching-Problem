#include "ILP_gurobi.h"
#include <string>
#include <iostream>
#include <cmath>
#include <limits>


std::vector<int> solveWithGurobi(const ProblemData& data) {
    std::vector<int> route;
    try {

        GRBEnv env(true);
        env.set("LogFile", "gurobi_tsp.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TSP_PrizeCollecting");

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
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                double objCoeff = (i == j) ? 0.0 : Prize[j];  
                x[i][j] = model.addVar(0.0, 1.0, objCoeff, GRB_BINARY,
                                       "x_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // 2) is_end[i]
        for(int i = 0; i < N; i++){
            is_end[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     "is_end_" + std::to_string(i));
        }

        // 3.3) u[i] 
        for(int i = 0; i < N; i++){
            if(i == start)
                continue;
            u[i] = model.addVar(2.0, (double)N, 0.0, GRB_INTEGER,
                                "u_" + std::to_string(i));
        }

        // 4) objective sense MAXIMIZE
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        // constraints
        // A) One_End_Node: sum_i is_end[i] = 1
        {
            GRBLinExpr lhs = 0.0;
            for(int i = 0; i < N; i++){
                lhs += is_end[i];
            }
            model.addConstr(lhs == 1, "One_End_Node");
        }

        // B) At_Most_One_Incoming: sum_{j != i} x[j,i] <= 1
        for(int i = 0; i < N; i++){
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++){
                if(j != i) lhs += x[j][i];
            }
            model.addConstr(lhs <= 1.0, "AtMostOneIncoming_" + std::to_string(i));
        }

        // C) At_Most_One_Outgoing: sum_{j != i} x[i,j] <= 1
        for(int i = 0; i < N; i++){
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++){
                if(j != i) lhs += x[i][j];
            }
            model.addConstr(lhs <= 1.0, "AtMostOneOutgoing_" + std::to_string(i));
        }

        // D) Flow_Balance (i != start): sum_j x[i,j] - sum_j x[j,i] + is_end[i] = 0
        for(int i = 0; i < N; i++){
            if(i == start) continue; 
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++){
                if(j != i) lhs += x[i][j];
            }
            for(int j = 0; j < N; j++){
                if(j != i) lhs -= x[j][i];
            }

            lhs += is_end[i];
            model.addConstr(lhs == 0.0, "FlowBalance_" + std::to_string(i));
        }

        // E) Start_Flow: sum_{j != start} x[start,j] - sum_{j != start} x[j,start] = 1
        {
            GRBLinExpr lhs = 0.0;
            for(int j = 0; j < N; j++){
                if(j != start) lhs += x[start][j];
            }
            for(int j = 0; j < N; j++){
                if(j != start) lhs -= x[j][start];
            }
            model.addConstr(lhs == 1.0, "StartFlow");
        }

        // F) Budget_Constraint: sum_{i != j} Cost[i][j]*x[i][j] <= Budget
        {
            GRBLinExpr lhs = 0.0;
            for(int i = 0; i < N; i++){
                for(int j = 0; j < N; j++){
                    if(i != j){
                        lhs += Cost[i][j]*x[i][j];
                    }
                }
            }
            model.addConstr(lhs <= Budget, "BudgetConstraint");
        }

        // G) Subtour_Elimination (MTZ):
        //    for i != j, i != start, j != start:
        //    u[i] - u[j] + N*x[i][j] <= N - 1
        for(int i = 0; i < N; i++){
            if(i == start) continue;
            for(int j = 0; j < N; j++){
                if(j == i || j == start) continue;
                GRBLinExpr lhs = 0.0;
                lhs += u[i];
                lhs -= u[j];
                lhs += (double)N * x[i][j];
                model.addConstr(lhs <= (double)N - 1.0,
                                "MTZ_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if(status == GRB_OPTIMAL) {
            double objVal = model.get(GRB_DoubleAttr_ObjVal);
            std::cout << "Optimal objective: " << objVal << "\n";
            
            // Build the route
            route.push_back(start);
            int currentNode = start;
            while(true) {
                bool foundNext = false;
                for(int nextNode = 0; nextNode < N; nextNode++){
                    if(nextNode == currentNode) continue;
                    double val = x[currentNode][nextNode].get(GRB_DoubleAttr_X);
                    if(val > 0.5){
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

            // // Print the route
            // std::cout << "Route found: ";
            // for(size_t i=0; i<route.size(); i++){
            //     std::cout << route[i];
            //     if(i < route.size() - 1){
            //         std::cout << " -> ";
            //     }
            // }
            // std::cout << "\n";


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


std::vector<int> gurobiSolve(const std::vector<std::vector<double>> &Cost, const std::vector<double> &Prize,
                int start, double Budget) {
    int N = Prize.size();
    ProblemData data(N, start, Budget, Cost, Prize);
    std::vector<int> path = solveWithGurobi(data);
    return path;
}
