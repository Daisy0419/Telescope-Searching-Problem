#include "TestFunctions.h"
#include "ReadData.h"
#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "ILP_gurobi.h"
#include "GCP.h"
#include "ILP_cplex.h"


#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem> 
#include <numeric>

//print path to stdout
double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> rank, std::vector<int> best_path, double padding) {

    if (best_path.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }

    if (rank.empty()) {
        rank.resize(best_path.size());
        std::iota(rank.begin(), rank.end(), 1); 
    }
    
    double sum_probability = 0.0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile << " "; 
        // std::cout << rank[tile] << " "; //ranks in tiling file
        sum_probability += probability[tile]; 
    }
    std::cout << std::endl;

    double total_cost = -padding;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        total_cost += costs[best_path[i]][best_path[i + 1]]; 
    }
    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; 
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}


void save_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, double slew_rate, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> rank,
                  std::vector<int> best_path, 
                  double elapsed_time, double padding) {

    bool file_exists = std::filesystem::exists(filename); 

    std::ofstream outfile(filename, std::ios::app);  

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    //Write header
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,SlewRate,NumTiles,SumProb,TotalCost,TimeSec,Path\n";
    }

    double sum_probability = 0.0;
    double total_cost = -padding;

    if (!best_path.empty()) {
        for (size_t i = 0; i < best_path.size(); i++) {
            sum_probability += probability[best_path[i]];
            if (i < best_path.size() - 1) {
                total_cost += costs[best_path[i]][best_path[i + 1]];
            }
        }
    }

    if (rank.empty()) {
        rank.resize(best_path.size());
        std::iota(rank.begin(), rank.end(), 1); 
    }

    // Write result to CSV
    outfile << method << ","
            << data << ","
            << budget << ","
            << slew_rate << ","
            << best_path.size() << ","
            << sum_probability << ","
            << total_cost << ","
            << elapsed_time << ",";

    // Write path
    for (size_t i = 0; i < best_path.size(); i++) {
        outfile << rank[best_path[i]];
        if (i < best_path.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << "\n";

    outfile.close();

}

// multi deadline
void save_result2(const std::string& filename, 
                  const std::string& method, 
                  const std::vector<double>& budgets,
                  const std::vector<double>& dwell_time, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> rank,
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
        outfile << "Method,NumTiles,TimeSec,Budgets,Costs,FOMs,Path\n";
    }

    double sum_probability = 0.0;
    double total_cost = 0.0;
    std::vector<double> foms;
    std::vector<double> accu_costs;
    int b = 0;

    if (!best_path.empty()) {
        for (size_t i = 1; i < best_path.size(); i++) {
            total_cost += costs[best_path[i]][best_path[i - 1]];

            if(b < budgets.size() && total_cost + 0.5*dwell_time[best_path[i]] > budgets[b]){
                foms.push_back(sum_probability);
                sum_probability = 0.0;
                double accumulated_cost=(total_cost - costs[best_path[i]][best_path[i - 1]] + 0.5*dwell_time[best_path[i-1]]);
                accu_costs.push_back(accumulated_cost);
                b++;
            }

            sum_probability += probability[best_path[i]];
        }
    }

    if (rank.empty()) {
        rank.resize(best_path.size());
        std::iota(rank.begin(), rank.end(), 1); 
    }

    // Write result to CSV
    outfile << method << ","
            << best_path.size() << ","
            << elapsed_time << ",";


    // Write budgets
    for (size_t i = 0; i < budgets.size(); i++) {
        outfile << budgets[i];
        if (i < budgets.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << ",";

    // Write Costs
    for (size_t i = 0; i < accu_costs.size(); i++) {
        outfile << accu_costs[i];
        if (i < accu_costs.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << ",";

    // Write foms
    for (size_t i = 0; i < foms.size(); i++) {
        outfile << foms[i];
        if (i < foms.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << ",";

    // Write path
    for (size_t i = 0; i < best_path.size(); i++) {
        outfile << rank[best_path[i]];
        if (i < best_path.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << "\n";

    outfile.close();

}


void test_algorithms_small_instances (std::string file, std::string out_file, 
                    double budget, double slew_rate, bool is_deepslow) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = 600;
    double accu_thr = 0.001;
    int init_pos_idx = 0;

    auto [start_idx, end_idx, padding] = buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                                                                slew_rate, is_deepslow, init_pos_idx);
         
    std::cout << "padding: " << padding << "\n";
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    
    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prizeGreedyPathTwoFixed(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Greedy", file, budget, slew_rate,costs, probability, ranks, greedy_path, elapsed_seconds.count(), padding);
    
    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int>  genetic_path = genetic_optimization_st(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Genetic", file, budget, slew_rate,costs, probability, ranks, genetic_path, elapsed_seconds.count(), padding);
    
    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs, probability);
    std::vector<int> ac_path = ac.ant_colony_optimization(budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, ac_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "AntColony", file, budget, slew_rate,costs, probability, ranks, ac_path, elapsed_seconds.count(), padding);

    // GCP
    std::cout << "*********running GCP*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_pathH = GCP(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_pathH, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Hoogeveen", file, budget, slew_rate,costs, probability, ranks, mst_pathH, elapsed_seconds.count(), padding);

  // gurobi with time limit
    std::cout << "*********gurobi solution with time limit*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> ilp_path = gurobiSolveST(costs, probability, start_idx, end_idx, budget+padding, accu_thr, time_limit, mst_pathH);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, ilp_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Gurobi", file, budget, slew_rate,costs, probability, ranks, ilp_path, elapsed_seconds.count(), padding);
    
}

void test_algorithms_large_instances (std::string file, std::string out_file, 
                    double budget, double slew_rate, bool is_deepslow) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = 1200;
    double accu_thr = 0.001;
    int init_pos_idx = 0;

    auto [start_idx, end_idx, padding] = buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                                                                slew_rate, is_deepslow, init_pos_idx);
         
    std::cout << "padding: " << padding << "\n";
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    
    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prizeGreedyPathTwoFixed(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Greedy", file, budget, slew_rate,costs, probability, ranks, greedy_path, elapsed_seconds.count(), padding);
    
    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int>  genetic_path = genetic_optimization_st(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Genetic", file, budget, slew_rate,costs, probability, ranks, genetic_path, elapsed_seconds.count(), padding);
    
    // GCP
    std::cout << "*********running GCP*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_pathH = GCP(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_pathH, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Hoogeveen", file, budget, slew_rate,costs, probability, ranks, mst_pathH, elapsed_seconds.count(), padding);

}

void test_algorithms_small_wcet(std::string file, std::string out_file, double budget, 
                            double budget_greedy, double budget_genetic, double budget_gcp,
                            double slew_rate, bool is_deepslow) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = -1;
    double accu_thr = 0.001;
    int init_pos_idx = 0;

    auto [start_idx, end_idx, padding] = buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                                                                slew_rate, is_deepslow, init_pos_idx);
         
    std::cout << "padding: " << padding << "\n";
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    
    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prizeGreedyPathTwoFixed(costs, probability, budget_greedy+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Greedy", file, budget, slew_rate, costs, probability, ranks, greedy_path, elapsed_seconds.count(), padding);
    
    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int>  genetic_path = genetic_optimization_st(costs, probability, budget_genetic+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Genetic", file, budget, slew_rate, costs, probability, ranks, genetic_path, elapsed_seconds.count(), padding);

    // mst Hoogeveen
    std::cout << "*********running GCP*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_pathH = GCP(costs, probability, budget_gcp+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_pathH, padding);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result(out_file, "Hoogeveen", file, budget, slew_rate, costs, probability, ranks, mst_pathH, elapsed_seconds.count(), padding);

}


void run_multi_deadlines(std::string file, std::string out_file, std::string out_file2, 
                                std::vector<double> budgets, 
                                double slew_rate, bool is_deepslow) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = 600;
    double accu_thr = -1;
    int init_pos_idx = 0;

    auto [start_idx, end_idx, padding] = buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                                                                slew_rate, is_deepslow, init_pos_idx);
         
    std::cout << "padding: " << padding << "\n";

    double budget1 = budgets[0], budget2 = budgets[1], budget3 = budgets[2];

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    double full_budget = budgets.back();

    for(int i = 0; i < budgets.size(); ++i) {
        double budget1 = budgets[i];
        double budget2 = full_budget;

        std::cout << "*********running start-Dn path*********" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::vector<int> d1_path = GCP(costs, probability, budget1+padding, start_idx, end_idx);

        std::cout << "*********running Dn-end path*********" << std::endl;
        std::vector<std::vector<double>> costs_remaining;
        std::vector<double> prizes_remaining;
        std::vector<int> ranks_mapping;
        std::vector<int> nodes_to_remove = d1_path;
        nodes_to_remove.pop_back();
        int next_start_idx = nodes_to_remove.back();

        auto [new_start_idx, new_end_idx] = removeNodesOrienteering(
            probability, costs, dwell_times, nodes_to_remove,
            costs_remaining, prizes_remaining, ranks_mapping,
            next_start_idx, end_idx
        );

        std::vector<int> d1_d2_path_ = GCP(
            costs_remaining, prizes_remaining,
            budget2 - budget1 + padding, new_start_idx, new_end_idx
        );

        std::vector<int> d1_d2_path = nodes_to_remove;
        for (int i = 1; i < d1_d2_path_.size(); ++i) {
            d1_d2_path.push_back(ranks_mapping[d1_d2_path_[i]]);
        }

        end = std::chrono::high_resolution_clock::now();
        print_path(costs, probability, ranks, d1_d2_path, padding);
        auto elapsed_seconds = end - start;

        std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

        std::string path_name = "D" + std::to_string(i);
        save_result(out_file, path_name, file, budget2 - budget1, slew_rate,
                    costs, probability, ranks, d1_d2_path, elapsed_seconds.count(), padding);

        save_result2(out_file2, path_name, budgets, dwell_times, costs, probability, ranks, 
                    d1_d2_path,  elapsed_seconds.count());
      
    }

}
