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
    std::string out_file  = "/home/research/w.yanwang/Telescope-Searching-Problem/results/RTSS_0518/out.csv";
    std::string out_file2 = "/home/research/w.yanwang/Telescope-Searching-Problem/results/RTSS_0518/out2.csv";
    std::string file      = "/home/research/w.yanwang/Telescope-Searching-Problem/data_RTSS_0518/large/filtered_GW191127_050227_7dt.csv";
    
    double slew_rate   = 50;
    bool is_deepslow   = false;
    int run_case       = 1;

    // Parse arguments
    if (argc > 1) file = std::string(argv[1]);
    if (argc > 3) {
        budget     = std::stod(argv[2]);
        run_case = std::stoi(argv[3]);
    }

    std::cout << "Input Parameters:\n";
    std::cout << "  file       = " << file       << "\n";
    std::cout << "  run_case   = " << run_case   << "\n\n";

    switch (run_case) {
        case 1:
            std::vector<double> budgets = {10, 20, 30, 40, 50, 60, 70, 80, 90};
            for(double budget in budgets)
                test_algorithms_small_instances(file, out_file, budget, slew_rate, is_deepslow);
            break;
        case 2:
            std::vector<double> budgets = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200};
            for(double budget in budgets)
                test_algorithms_large_instances(file, out_file, budget, slew_rate, is_deepslow);
            break;
        case 3: {
            std::vector<double> budgets = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200};
            run_multi_deadlines(file, out_file, out_file2, budgets, slew_rate, dwell_time, is_deepslow);
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


