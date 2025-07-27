#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "ILP_gurobi.h"
#include "GCP.h"
#include "ILP_cplex.h"
#include "TestFunctions.h"


int main(int argc, char** argv) {
    // Default parameters
    std::string out_file  = "../results/recomputed_results/out.csv";
    std::string out_file2 = "../results/recomputed_results/out2.csv";
    std::string file      = "../data/small/filtered_GW191105_143521_7dt_separate.csv";
    
    double slew_rate   = 50;
    int run_case       = 1;
    bool is_deepslow   = false;
    double budget      = 50;

    // Parse arguments
    if (argc > 1) file = std::string(argv[1]);
    if (argc > 2) {
        run_case = std::stoi(argv[2]);
    }
    if (argc > 3) {
        budget = std::stod(argv[3]);
    }

    double budget_greedy=1, budget_genetic=1, budget_gcp=1;
    if (argc > 6) {
        budget_greedy = std::stod(argv[4]);
        budget_genetic = std::stod(argv[5]);
        budget_gcp = std::stod(argv[6]);
    }

    std::cout << "Input Parameters:\n";
    std::cout << "  file       = " << file       << "\n";
    std::cout << "  run_case   = " << run_case   << "\n\n";

    switch (run_case) {
        case 1:
            test_algorithms_small_instances(file, out_file, budget, slew_rate, is_deepslow);
            break;
        case 2:
            test_algorithms_large_instances(file, out_file, budget, slew_rate, is_deepslow);
            break;
        case 3:
            test_algorithms_small_wcet(file, out_file, budget, budget_greedy, budget_genetic, 
                                        budget_gcp, slew_rate, is_deepslow);
            break;
        case 4: {
            std::vector<double> budgets = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200};
            run_multi_deadlines(file, out_file, out_file2, budgets, slew_rate, is_deepslow);
            break;
        }
        default:
            std::cerr << "Invalid run_case. Use:\n"
                      << "  1 = small instances\n"
                      << "  2 = large instances\n"
                      << "  3 = multi-deadline instances\n";
    }
    return 0;
}


