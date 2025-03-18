#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "Christofides.h"
#include "kMST.h"
#include "ReadData.h"
#include "BranchBound.h"
// #include "ILP_scip.h"
#include "ILP_gurobi.h"

#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem> 

double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                std::vector<int> best_path, double dwell_time) {

    if (best_path.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }
    
    double sum_probability = 0.0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile + 1 << " "; //convert 1-indexing
        sum_probability += probability[tile]; // Accumulate probability
    }
    std::cout << std::endl;

    double total_cost = dwell_time;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        total_cost += costs[best_path[i]][best_path[i + 1]]; // Compute total cost
    }

    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; // Correct value
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}


void write_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, double slew_rate, double dwell_time, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> best_path, 
                  double elapsed_time) {

    bool file_exists = std::filesystem::exists(filename); 

    std::ofstream outfile(filename, std::ios::app);  

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    //Write header
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,SlewRate,DwellTime,NumTiles,SumProb,TotalCost,TimeSec,Path\n";
    }

    double sum_probability = 0.0;
    double total_cost = dwell_time;

    if (!best_path.empty()) {
        for (size_t i = 0; i < best_path.size(); i++) {
            sum_probability += probability[best_path[i]];
            if (i < best_path.size() - 1) {
                total_cost += costs[best_path[i]][best_path[i + 1]];
            }
        }
    }

    // Write result to CSV
    outfile << method << ","
            << data << ","
            << budget << ","
            << slew_rate << ","
            << dwell_time << ","
            << best_path.size() << ","
            << sum_probability << ","
            << total_cost << ","
            << elapsed_time << ",\"";

    // Write path (space-separated inside quotes)
    for (size_t i = 0; i < best_path.size(); i++) {
        outfile << best_path[i] + 1;
        if (i < best_path.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << "\"\n";

    outfile.close();

    // // Print to console for debugging (optional)
    // std::cout << "[RESULT] Method: " << method
    //           << " | Dataset: " << file
    //           << " | Budget: " << budget
    //           << " | SlewRate: " << slew_rate
    //           << " | DwellTime: " << dwell_time
    //           << " | Path: ";
    // for (int tile : best_path) std::cout << tile << " ";
    // std::cout << "| NumTiles: " << best_path.size()
    //           << " | SumProb: "  << sum_probability
    //           << " | TotalCost: " << total_cost
    //           << " | TimeSec: "   << elapsed_time
    //           << std::endl;
}


void test_all_algorithms (std::string file, std::string out_file, double budget, 
                            double slew_rate, double dwell_time,  bool is_deepslow=false) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;

    if(is_deepslow) {
        read_data_deep_slow(file, costs, probability, slew_rate, dwell_time);
    } else {
        read_data(file, costs, probability, slew_rate, dwell_time);
    }
    
    // read_data(file, costs, probability, slew_rate, dwell_time);    
        
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs, probability);
    std::vector<int> ac_path = ac.ant_colony_optimization(budget-dwell_time, 0);
    end = std::chrono::high_resolution_clock::now();
    double best_prize = print_path(costs, probability, ac_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "AntColony", file, budget, slew_rate, dwell_time, costs, probability, ac_path, elapsed_seconds.count());

    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int> genetic_path = genetic_optimization(costs, probability, budget-dwell_time, 0);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, genetic_path, elapsed_seconds.count());


    // //greedy
    // std::cout << "*********running greedy*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // std::vector<int> greedy_path = prize_greedy_path(costs, probability, 0, budget-dwell_time);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, greedy_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, greedy_path, elapsed_seconds.count());


    // // branch and bound
    // std::cout << "*********running branch and bound*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // std::vector<int> true_path = run_branch_bound(costs, probability, budget-dwell_time, 0, ac_path, best_prize);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, true_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    // //gurobi
    // std::cout << "*********running gurobi solver*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time);
    // // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time, genetic_path);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, ilp_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "Gurobi", file, budget, slew_rate, dwell_time, costs, probability, ilp_path, elapsed_seconds.count());
    
    // //kmst
    // std::cout << "*********running kmst*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // findBestTour(costs, budget);
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

}


void test_combined_algorithms (std::string file, std::string out_file, double budget, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<std::vector<double>> region_costs;
    std::vector<double> region_probability;
    std::vector<double> ranks;
    double region_sum_probability = 0.99;
    double time_limit = 30;
    double accu_thr = 0.001;

    if(is_deepslow) {
        // read_data_deep_slow(file, costs, probability, slew_rate, dwell_time);
        read_data_deep_slow_filtered(file, costs, probability, region_costs, region_probability, ranks, 
                slew_rate, dwell_time, region_sum_probability);
    } else {
        // read_data(file, costs, probability, slew_rate, dwell_time);
        read_data_filtered(file, costs, probability, region_costs, region_probability, ranks, 
                        slew_rate, dwell_time, region_sum_probability);
    }

    std::cout << "region_probability.size(): " << region_probability.size() << std::endl;
      
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    // gurobi with time limit
    std::cout << "*********gurobi solution with time/accuracy restrain*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time, accu_thr, time_limit, {});
    // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time, genetic_path);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ilp_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Gurobi+Restrain", file, budget, slew_rate, dwell_time, costs, probability, ilp_path, elapsed_seconds.count());
    
    //regional solution with gurobi
    std::cout << "*********get regional solution with gurobi solver*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    ilp_path = gurobiSolve(region_costs, region_probability, 0, budget-dwell_time, accu_thr, time_limit, {});
    // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time, genetic_path);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ilp_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Gurobi+Region", file, budget, slew_rate, dwell_time, costs, probability, ilp_path, elapsed_seconds.count());
    

    //combined solution with gurobi + greedy
    std::cout << "*********gurobi + greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    double region_budget = 0.9 * budget;
    ilp_path = gurobiSolve(region_costs, region_probability, 0, region_budget-dwell_time, accu_thr, time_limit, {});
    std::vector<int> greedy_path = partial_prize_greedy_path(costs, probability, ilp_path, ilp_path.back(), budget-region_budget);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Gurobi+Greedy", file, budget, slew_rate, dwell_time, costs, probability, greedy_path, elapsed_seconds.count());


    //combined solution with gurobi + genetric
    std::cout << "*********gurobi + genetric*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    region_budget = budget;
    ilp_path = gurobiSolve(region_costs, region_probability, 0, region_budget-dwell_time, accu_thr, time_limit, {});
    std::vector<std::vector<int>> gene_initial_paths;
    for(size_t i = 0; i < 10; ++i) {
        gene_initial_paths.push_back(ilp_path);
    }
    std::vector<int> genetic_path = genetic_optimization(costs, probability, budget-dwell_time, 0, gene_initial_paths);

    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Gurobi+Genetic", file, budget, slew_rate, dwell_time, costs, probability, genetic_path, elapsed_seconds.count());

    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    greedy_path = prize_greedy_path(costs, probability, 0, budget-dwell_time);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, greedy_path, elapsed_seconds.count());


    // // ant colony
    // std::cout << "*********running ant colony*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // AntColony ac = AntColony(costs, probability);
    // std::vector<int> ac_path = ac.ant_colony_optimization(budget-dwell_time, 0);
    // end = std::chrono::high_resolution_clock::now();
    // double best_prize = print_path(costs, probability, ac_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "AntColony", file, budget, slew_rate, dwell_time, costs, probability, ac_path, elapsed_seconds.count());

    // genetic + greedy
    std::cout << "*********running genetic + greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int>> greedy_initial_paths;
    for(size_t i = 0; i < 10; ++i) {
        greedy_initial_paths.push_back(greedy_path);
    }
    genetic_path = genetic_optimization(costs, probability, budget-dwell_time, 0, greedy_initial_paths);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy+Genetic", file, budget, slew_rate, dwell_time, costs, probability, genetic_path, elapsed_seconds.count());


    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    genetic_path = genetic_optimization(costs, probability, budget-dwell_time, 0);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, genetic_path, elapsed_seconds.count());

}

int main (int argc, char** argv) {
    // default paramaters
    std::string out_file = "../results/out.csv";
    // std::string file = "../data/7dt_combined.csv";
    // std::string file = "../data/7dt_separate.csv";
    std::string file = "../data/deep_slow.csv";
    
    double budget = 100;
    double slew_rate = 1;
    double dwell_time = 0.0;
    bool is_deepslow = true;

    if (argc > 1) {
        file = std::string(argv[1]);
    }
    if (argc > 4) {
        budget = std::stod(argv[2]);
        slew_rate = std::stod(argv[3]);
        dwell_time = std::stod(argv[4]);
        is_deepslow = std::stoi(argv[5]);
    }

    std::cout << "Input Parameters:\n";
    std::cout << "  file       = " << file       << "\n";
    std::cout << "  budget     = " << budget     << "\n";
    std::cout << "  slew_rate  = " << slew_rate  << "\n";
    std::cout << "  dwell_time = " << dwell_time << "\n\n";


    test_combined_algorithms (file, out_file, budget, slew_rate, dwell_time, is_deepslow);
    // test_all_algorithms (file, out_file, budget, slew_rate, dwell_time, is_deepslow);
    return 0;
}


