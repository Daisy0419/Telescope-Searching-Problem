// #include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "ILP_gurobi.h"
#include "GCP.h"
#include "TestFunctions.h"

#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>


static std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

template <class F>
static std::vector<int> run_and_report(const char* tag, F&& f,
                                       const std::vector<std::vector<double>>& costs,
                                       const std::vector<double>& probability,
                                       const std::vector<int>& ranks,
                                       double padding)
{
    using Clock = std::chrono::high_resolution_clock;
    std::cout << "********* running " << tag << " *********\n";
    auto t0 = Clock::now();
    std::vector<int> path = f();
    auto t1 = Clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;

    print_path(costs, probability, ranks, path, padding);
    std::cout << "running time (wallclock): " << elapsed.count() << " seconds\n";
    return path;
}

int main(int argc, char** argv) {
    // Defaults
    std::string out_file  = "../results/recomputed_results/out.csv";
    std::string out_file2 = "../results/recomputed_results/out2.csv";
    std::string file      = "../data/small/filtered_GW191105_143521_7dt_separate.csv";

    double slew_rate   = 50.0;
    bool is_deepslow   = false;
    double budget      = 50.0;
    std::string alg    = "greedy";

    // Parse: file budget alg [slew_rate] [is_deepslow]
    if (argc < 4) {
        std::cerr <<
            "Usage:\n"
            "  " << argv[0] << " <file> <budget> <alg> [slew_rate=50]\n"
            "Where <alg> in: greedy | genetic | gcp | ilp \n";
        return 1;
    } else {
        file   = std::string(argv[1]);
        budget = std::stod(argv[2]);
        alg    = std::string(argv[3]);
        if (argc > 4) slew_rate   = std::stod(argv[4]);
    }

    std::string alg_lc = lower(alg);


    // Build graph
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;

    double time_limit = 600.0;
    double accu_thr   = 0.001;
    int init_pos_idx  = 0;

    auto [start_idx, end_idx, padding] =
        buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                               slew_rate, is_deepslow, init_pos_idx);

    std::cout << "Input Parameters:\n";
    std::cout << "  file       = " << file       << "\n";
    std::cout << "  budget     = " << budget     << "\n";
    std::cout << "  alg        = " << alg        << "\n";
    std::cout << "  slew_rate  = " << slew_rate  << "\n";

    const double eff_budget = budget + padding;

    if (alg_lc == "greedy") {
        (void)run_and_report("Greedy", [&]{
            return prizeGreedyPathTwoFixed(costs, probability, eff_budget, start_idx, end_idx);
        }, costs, probability, ranks, padding);
    }
    else if (alg_lc == "genetic") {
        (void)run_and_report("Genetic", [&]{
            return genetic_optimization_st(costs, probability, eff_budget, start_idx, end_idx);
        }, costs, probability, ranks, padding);
    }
    else if (alg_lc == "gcp") {
        (void)run_and_report("GCP", [&]{
            return GCP(costs, probability, eff_budget, start_idx, end_idx);
        }, costs, probability, ranks, padding);
    }
    else if (alg_lc == "ilp" || alg_lc == "gurobi") {
        // Optional warm start via GCP
        // std::vector<int> warm = GCP(costs, probability, eff_budget, start_idx, end_idx);
        (void)run_and_report("Gurobi", [&]{
            return gurobiSolveST(costs, probability, start_idx, end_idx,
                                 eff_budget, accu_thr, time_limit, {});
        }, costs, probability, ranks, padding);
    }
    else {
        std::cerr << "Please input a valid algorithm: greedy, genetic, gcp, ilp, or aco\n";
        return 2;
    }

    return 0;
}
