#include "ILP_scip.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"


//constraints and variables
SCIP_RETCODE tspProblem::scip_create(SCIP* scip) {
    SCIPdebugMessage("tspProblem::scip_create() called.\n");

    int N = data->N;
    int startIndex = data->startIndex;
    double Budget = data->Budget;
    const auto& Cost  = data->Cost;
    const auto& Prize = data->Prize;

    // 1) x[i][j], BINARY
    x.resize(N);
    for(int i=0; i<N; i++){
        x[i].resize(N, nullptr);
        for(int j=0; j<N; j++){
            double objCoeff = (i == j) ? 0.0 : Prize[j];

            std::string varName = "x_" + std::to_string(i) + "_" + std::to_string(j);
            SCIP_VAR* var = nullptr;
            SCIP_CALL( SCIPcreateVarBasic(
                scip,
                &var,
                varName.c_str(),
                0.0, 1.0,   
                objCoeff,  
                SCIP_VARTYPE_BINARY
            ) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            x[i][j] = var;
        }
    }

    // 2) is_end[i], BINARY
    is_end.resize(N, nullptr);
    for(int i=0; i<N; i++){
        std::string vname = "is_end_" + std::to_string(i);
        SCIP_VAR* var = nullptr;
        SCIP_CALL( SCIPcreateVarBasic(
            scip,
            &var,
            vname.c_str(),
            0.0, 1.0,
            0.0,  
            SCIP_VARTYPE_BINARY
        ));
        SCIP_CALL( SCIPaddVar(scip, var) );
        is_end[i] = var;
    }

    // 3) u[i]
    u.resize(N, nullptr);
    for(int i=1; i<N; i++){
        std::string vname = "u_" + std::to_string(i);
        SCIP_VAR* var = nullptr;
        SCIP_CALL( SCIPcreateVarBasic(
            scip,
            &var,
            vname.c_str(),
            2.0, (SCIP_Real)N, 
            0.0,
            SCIP_VARTYPE_INTEGER 
        ));
        SCIP_CALL( SCIPaddVar(scip, var) );
        u[i] = var;
    }

    // 4) Constraints
    // A) One_End_Node: sum_i is_end[i] = 1
    {
        SCIP_CONS* cons = nullptr;
        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, "One_End_Node",
            0, nullptr, nullptr,
            1.0, 1.0  // lhs=1, rhs=1
        ));
        for(int i=0; i<N; i++){
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, is_end[i], 1.0) );
        }
        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // B) At_Most_One_Incoming: for each i, sum_{j != i} x[j,i] <= 1
    for(int i=0; i<N; i++){
        std::string cname = "At_Most_One_Incoming_" + std::to_string(i);
        SCIP_CONS* cons = nullptr;
        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, cname.c_str(),
            0, nullptr, nullptr,
            0, 1.0
        ));
        for(int j=0; j<N; j++){
            if(j != i){
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[j][i], 1.0) );
            }
        }
        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // C) At_Most_One_Outgoing: for each i, sum_{j != i} x[i,j] <= 1
    for(int i=0; i<N; i++){
        std::string cname = "At_Most_One_Outgoing_" + std::to_string(i);
        SCIP_CONS* cons = nullptr;
        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, cname.c_str(),
            0, nullptr, nullptr,
            0, 1.0
        ));
        for(int j=0; j<N; j++){
            if(j != i){
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i][j], 1.0) );
            }
        }
        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // D) Connectivity (i != start): sum_j x[i,j] - sum_j x[j,i] = -is_end[i]
    for(int i=0; i<N; i++){
        if(i == startIndex) 
            continue;
        std::string cname = "Flow_Balance_" + std::to_string(i);
        SCIP_CONS* cons = nullptr;

        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, cname.c_str(),
            0, nullptr, nullptr,
            0.0, 0.0
        ));
        // sum_j x[i,j]
        for(int j=0; j<N; j++){
            if(j != i){
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i][j],  -1.0) );
            }
        }
        // - sum_j x[j,i]
        for(int j=0; j<N; j++){
            if(j != i){
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[j][i], 1.0) );
            }
        }
        // - is_end[i]
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, is_end[i], -1.0) );

        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // E) Start_Flow: sum_{j != start} x[start,j] - sum_{j != start} x[j,start] = 1
    {
        std::string cname = "Start_Flow";
        SCIP_CONS* cons = nullptr;
        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, cname.c_str(),
            0, nullptr, nullptr,
            1.0, 1.0
        ));
        for(int j=0; j<N; j++){
            if(j != startIndex){
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[startIndex][j],  1.0) );
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[j][startIndex], -1.0) );
            }
        }
        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // F) Budget_Constraint: sum_{i != j} Cost[i][j]*x[i][j] <= Budget
    {
        SCIP_CONS* cons = nullptr;
        SCIP_CALL( SCIPcreateConsBasicLinear(
            scip, &cons, "Budget_Constraint",
            0, nullptr, nullptr,
            0, Budget
        ));
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                if(i != j){
                    SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i][j], Cost[i][j]) );
                }
            }
        }
        SCIP_CALL( SCIPaddCons(scip, cons) );
    }

    // G) Subtour_Elimination 
    //    For i != j, i != start, j != start:
    //    u[i] - u[j] + N * x[i][j] <= N - 1

    for (int i = 0; i < N; i++) {
        if (i == startIndex) continue;
        for (int j = 0; j < N; j++) {
            if (j == i || j == startIndex) continue;

            std::string cname = "MTZ_" + std::to_string(i) + "_" + std::to_string(j);
            SCIP_CONS* cons = nullptr;

            std::vector<SCIP_VAR*> vars;
            std::vector<SCIP_Real> coefs;

            // u[i] * 1
            vars.push_back(u[i]);
            coefs.push_back(1.0);

            // u[j] * -1
            vars.push_back(u[j]);
            coefs.push_back(-1.0);

            // N * x[i][j]
            vars.push_back(x[i][j]);
            coefs.push_back((SCIP_Real)N);

            SCIP_VAR** vars_array = vars.data();
            SCIP_Real* coefs_array = coefs.data();

            SCIP_CALL(SCIPcreateConsBasicLinear(
                scip, &cons, cname.c_str(),
                vars.size(), vars_array, coefs_array,
                -(SCIP_Real)N + 1.0, (SCIP_Real)N - 1.0  // lhs , rhs
            ));

            SCIP_CALL(SCIPaddCons(scip, cons));
        }
    }

    return SCIP_OKAY;
}

int solveProblem(const std::vector<std::vector<double>> &Cost, const std::vector<double> &Prize,
                double Budget, int startIndex) {
    SCIP* scip = nullptr;
    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    int N = Prize.size();
    auto myData = std::make_shared<problemData>(N, startIndex, Budget, Cost, Prize);
    auto myProb = std::make_shared<tspProblem>(myData);

    for(int i = 0; i < N; ++i) {
        std::cout << (myData->Prize)[i] << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++ j) {
            std::cout << (myData->Cost)[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "budget: " << myData->Budget << std::endl;
    std::cout << "startIndex: " << myData->startIndex << std::endl;
    std::cout << "N: " << myData->N << std::endl;
    std::cout << "myProb.get(): " << myProb.get() << std::endl;
    

    SCIP_CALL(SCIPcreateObjProb(scip, "tsp", myProb.get(), false));
    SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);

    SCIP_CALL(SCIPsolve(scip));
    SCIP_SOL* sol = SCIPgetBestSol(scip);
    
    // if (sol != nullptr) {
    //     std::cout << "Optimal objective: " << SCIPgetSolOrigObj(scip, sol) << "\n";

    //     std::vector<int> route;
    //     double totalCost = 0.0;

    //     int currentNode = startIndex;
    //     route.push_back(currentNode);

    //     while (true) {
    //         bool foundNext = false;
    //         for (int nextNode = 0; nextNode < N; ++nextNode) {
    //             if (currentNode != nextNode && 
    //                 SCIPgetSolVal(scip, sol, myProb->x[currentNode][nextNode]) > 0.5) {

    //                 route.push_back(nextNode);
    //                 totalCost += Cost[currentNode][nextNode];
    //                 currentNode = nextNode;
    //                 foundNext = true;
    //                 break;
    //             }
    //         }

    //         if (!foundNext) {
    //             break;
    //         }
    //     }

    //     std::cout << "Found route: ";
    //     for (size_t i = 0; i < route.size(); ++i) {
    //         std::cout << route[i];
    //         if (i < route.size() - 1) {
    //             std::cout << " -> ";
    //         }
    //     }
    //     std::cout << "\n";

    //     std::cout << "Total cost of the route: " << totalCost << "\n";

    // } else {
    //     std::cout << "No solution found.\n";
    // }

    SCIP_CALL(SCIPfree(&scip));
    return 0;
}


