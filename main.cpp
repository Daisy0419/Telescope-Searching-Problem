#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "BranchBound.h"
#include "ILP_gurobi.h"
#include "SimulatedAnnealing.h"
#include "bpc_tsp.h"
#include "mst.h"

#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem> 
#include <numeric>

double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> rank, std::vector<int> best_path, double dwell_time) {

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
        std::cout << rank[tile] << " "; 
        sum_probability += probability[tile]; 
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    // std::cout << " segmemt costs: ";
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        // std::cout << best_path[i] << " " << best_path[i + 1]<< " " << costs[best_path[i]][best_path[i + 1]] << ", ";
        total_cost += costs[best_path[i]][best_path[i + 1]]; 
    }
    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; 
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}

void write_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, double slew_rate, double dwell_time, 
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
        outfile << "Method,Dataset,Budget,SlewRate,DwellTime,NumTiles,SumProb,TotalCost,TimeSec,Path\n";
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
            << dwell_time << ","
            << best_path.size() << ","
            << sum_probability << ","
            << total_cost << ","
            << elapsed_time << ",";

    // Write path
    // outfile << rank[init_pos_idx] << " ";
    for (size_t i = 0; i < best_path.size(); i++) {
        outfile << rank[best_path[i]];
        if (i < best_path.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << "\n";

    outfile.close();

}

void write_result2(const std::string& filename, 
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
    // outfile << rank[init_pos_idx] << " ";
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

void test_algorithms (std::string file, std::string out_file, double budget, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {

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
    std::vector<int> greedy_path = prizeGreedyPathTwoFixed(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, ranks, greedy_path, elapsed_seconds.count(), padding);

    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int>  genetic_path = genetic_optimization_st(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, ranks, genetic_path, elapsed_seconds.count(), padding);
    
    // // ant colony
    // std::cout << "*********running ant colony*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // AntColony ac = AntColony(costs, probability);
    // std::vector<int> ac_path = ac.ant_colony_optimization(budget+padding, start_idx, end_idx);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, ranks, ac_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "AntColony", file, budget, slew_rate, dwell_time, costs, probability, ranks, ac_path, elapsed_seconds.count(), padding);


    // mst Hoogeveen
    std::cout << "*********running Hoogeveen*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_pathH = bestHoogeveen(costs, probability, budget+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_pathH, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Hoogeveen", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_pathH, elapsed_seconds.count(), padding);


//     // mst Hoogeveen
//     std::cout << "*********running HoogeveenParallel*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> mst_pathP = bestHoogeveenParallel(costs, probability, budget+padding, start_idx, end_idx);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, mst_pathP, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "HoogeveenParallel", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_pathP, elapsed_seconds.count(), padding);

//     // greedy + genetic
//     std::cout << "*********running Hoogeveen+Genetic*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<std::vector<int>> initial_paths;
//     for(size_t i = 0; i < 20; ++i) {
//         initial_paths.push_back(mst_pathH);
//     }
//     std::vector<int> genetic_path2 = genetic_optimization_st(costs, probability, budget+padding, start_idx, end_idx, initial_paths);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, genetic_path2, dwell_time);
//     elapsed_seconds = elapsed_seconds + end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "Hoogeveen+Genetic", file, budget, slew_rate, dwell_time, costs, probability, ranks, genetic_path2, elapsed_seconds.count(), padding);

//   // gurobi with time limit
//     std::cout << "*********gurobi solution with time limit*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> ilp_path = gurobiSolveST(costs, probability, start_idx, end_idx, budget+padding, accu_thr, time_limit);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, ilp_path, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "Gurobi", file, budget, slew_rate, dwell_time, costs, probability, ranks, ilp_path, elapsed_seconds.count(), padding);

}


void test_algorithms2(std::string file, std::string out_file, double budget, 
                            double budget_greedy, double budget_genetic, double budget_gcp,
                            double slew_rate, double dwell_time, bool is_deepslow=false) {

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
    print_path(costs, probability, ranks, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, ranks, greedy_path, elapsed_seconds.count(), padding);

    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int>  genetic_path = genetic_optimization_st(costs, probability, budget_genetic+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, ranks, genetic_path, elapsed_seconds.count(), padding);

    // mst Hoogeveen
    std::cout << "*********running Hoogeveen*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_pathH = bestHoogeveen(costs, probability, budget_gcp+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_pathH, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Hoogeveen", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_pathH, elapsed_seconds.count(), padding);

}


void test_multi_deadlines(std::string file, std::string out_file, std::string out_file2, 
                                std::vector<double> budgets, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = 900;
    double accu_thr = -1;
    int init_pos_idx = 0;

    auto [start_idx, end_idx, padding] = buildGraphOrienteering(file, costs, probability, ranks, dwell_times,
                                                                slew_rate, is_deepslow, init_pos_idx);
         
    std::cout << "padding: " << padding << "\n";

    double budget1 = budgets[0], budget2 = budgets[1], budget3 = budgets[2];

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    std::cout << "*********running d1 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d1_path = bestHoogeveen(costs, probability, budget1+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d1_path, dwell_time);
    auto elapsed_seconds_d1 = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds_d1.count() << "seconds" << std::endl;
    write_result(out_file, "d1", file, budget1, slew_rate, dwell_time, costs, probability, ranks, d1_path, elapsed_seconds_d1.count(), padding);

    std::cout << "*********running d2 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d2_path = bestHoogeveen(costs, probability, budget2+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d2_path, dwell_time);
    auto elapsed_seconds_d2 = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds_d2.count() << "seconds" << std::endl;
    write_result(out_file, "d2", file, budget2, slew_rate, dwell_time, costs, probability, ranks, d2_path, elapsed_seconds_d2.count(), padding);

    std::cout << "*********running d3 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d3_path = bestHoogeveen(costs, probability, budget3+padding, start_idx, end_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d3_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d3", file, budget3, slew_rate, dwell_time, costs, probability, ranks, d3_path, elapsed_seconds.count(), padding);
    write_result2(out_file2, "d3", budgets, dwell_times, costs, probability, ranks, d3_path,  elapsed_seconds.count());
      
    std::vector<std::vector<double>> costs_remaining;
    std::vector<double> prizes_remaining;
    std::vector<int> ranks_mapping;
    std::cout << "*********running d1-d2 path*********" << std::endl;

    start = std::chrono::high_resolution_clock::now();

    std::vector<int> nodes_to_remove = d1_path;
    nodes_to_remove.pop_back();
    int next_start_idx = nodes_to_remove.back();

    auto [new_start_idx, new_end_idx] = removeNodesOrienteering(
        probability, costs, dwell_times, nodes_to_remove,
        costs_remaining, prizes_remaining, ranks_mapping,
        next_start_idx, end_idx
    );

    std::vector<int> d1_d2_path_ = bestHoogeveen(
        costs_remaining, prizes_remaining,
        budget2 - budget1 + padding, new_start_idx, new_end_idx
    );

    std::vector<int> d1_d2_path = nodes_to_remove;
    for (int i = 1; i < d1_d2_path_.size(); ++i) {
        d1_d2_path.push_back(ranks_mapping[d1_d2_path_[i]]);
    }

    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d1_d2_path, dwell_time);
    auto elapsed_seconds_d1_d2 = end - start + elapsed_seconds_d1;

    std::cout << "running time (wallclock): " << elapsed_seconds_d1_d2.count() << " seconds" << std::endl;
    write_result(out_file, "d1-d2", file, budget2 - budget1, slew_rate, dwell_time,
                costs, probability, ranks, d1_d2_path, elapsed_seconds_d1_d2.count(), padding);


    std::cout << "*********running d2-d3 path*********" << std::endl;
    {
        start = std::chrono::high_resolution_clock::now();

        std::vector<int> nodes_to_remove2 = d2_path;
        nodes_to_remove2.pop_back();
        int next_start_idx = nodes_to_remove2.back();

        auto [new_start_idx, new_end_idx] = removeNodesOrienteering(
            probability, costs, dwell_times, nodes_to_remove2,
            costs_remaining, prizes_remaining, ranks_mapping,
            next_start_idx, end_idx
        );

        std::vector<int> d2_d3_path_ = bestHoogeveen(
            costs_remaining, prizes_remaining,
            budget3 - budget2 + padding, new_start_idx, new_end_idx
        );

        std::vector<int> d2_d3_path = nodes_to_remove2;
        for (int i = 1; i < d2_d3_path_.size(); ++i) {
            d2_d3_path.push_back(ranks_mapping[d2_d3_path_[i]]);
        }

        end = std::chrono::high_resolution_clock::now();
        print_path(costs, probability, ranks, d2_d3_path, dwell_time);
        elapsed_seconds = end - start + elapsed_seconds_d2;

        std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;
        write_result(out_file, "d2_d3", file, budget3 - budget2, slew_rate, dwell_time,
                    costs, probability, ranks, d2_d3_path, elapsed_seconds.count(), padding);
        write_result2(out_file2, "d2_d3", budgets, dwell_times, costs, probability, ranks, d2_d3_path,  elapsed_seconds.count());
        
    }
    std::cout << "*********running d1-d2-d3 path*********" << std::endl;
    {
        start = std::chrono::high_resolution_clock::now();

        std::vector<int> nodes_to_remove3 = d1_d2_path;
        nodes_to_remove3.pop_back();
        int next_start_idx = nodes_to_remove3.back();

        auto [new_start_idx, new_end_idx] = removeNodesOrienteering(
            probability, costs, dwell_times, nodes_to_remove3,
            costs_remaining, prizes_remaining, ranks_mapping,
            next_start_idx, end_idx
        );

        std::vector<int> d1_d2_d3_path_ = bestHoogeveen(
            costs_remaining, prizes_remaining,
            budget3 - budget2 + padding, new_start_idx, new_end_idx
        );

        std::vector<int> d1_d2_d3_path = nodes_to_remove3;
        for (int i = 1; i < d1_d2_d3_path_.size(); ++i) {
            d1_d2_d3_path.push_back(ranks_mapping[d1_d2_d3_path_[i]]);
        }

        end = std::chrono::high_resolution_clock::now();
        print_path(costs, probability, ranks, d1_d2_d3_path, dwell_time);
        elapsed_seconds = end - start + elapsed_seconds_d1_d2;

        std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;
        write_result(out_file, "d1_d2_d3", file, budget3 - budget2, slew_rate, dwell_time,
                    costs, probability, ranks, d1_d2_d3_path, elapsed_seconds.count(), padding);
        write_result2(out_file2, "d1_d2_d3", budgets, dwell_times, costs, probability, ranks, d1_d2_d3_path,  elapsed_seconds.count());
      
    }

}



void run_multi_deadlines(std::string file, std::string out_file, std::string out_file2, 
                                std::vector<double> budgets, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
   
    double time_limit = 900;
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
        std::vector<int> d1_path = bestHoogeveen(costs, probability, budget1+padding, start_idx, end_idx);

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

        std::vector<int> d1_d2_path_ = bestHoogeveen(
            costs_remaining, prizes_remaining,
            budget2 - budget1 + padding, new_start_idx, new_end_idx
        );

        std::vector<int> d1_d2_path = nodes_to_remove;
        for (int i = 1; i < d1_d2_path_.size(); ++i) {
            d1_d2_path.push_back(ranks_mapping[d1_d2_path_[i]]);
        }

        end = std::chrono::high_resolution_clock::now();
        print_path(costs, probability, ranks, d1_d2_path, dwell_time);
        auto elapsed_seconds = end - start;

        std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

        std::string path_name = "D" + std::to_string(i);
        write_result(out_file, path_name, file, budget2 - budget1, slew_rate, dwell_time,
                    costs, probability, ranks, d1_d2_path, elapsed_seconds.count(), padding);

        write_result2(out_file2, path_name, budgets, dwell_times, costs, probability, ranks, 
                    d1_d2_path,  elapsed_seconds.count());
      
    }

}


int main (int argc, char** argv) {
    // default paramaters
    std::string out_file = "/home/research/w.yanwang/Telescope-Searching-Problem/results/RTSS_0518/out.csv";
    std::string out_file2 = "/home/research/w.yanwang/Telescope-Searching-Problem/results/RTSS_0518/out2.csv";
    std::string file = "../data_RTSS/large_cases/filtered_GW200105_162426.fits_7dt.csv";//59 64 70 83 91 104
    double budget = 1000;
    double slew_rate = 50;
    double dwell_time = 1;
    bool is_deepslow = false;
    
    if (argc > 1) {
        file = std::string(argv[1]);
    }
    if (argc > 2) {
        budget = std::stod(argv[2]);
        // slew_rate = std::stod(argv[3]);
        // dwell_time = std::stod(argv[4]);
        // is_deepslow = std::stoi(argv[5]);
    }
    double budget_greedy=1, budget_genetic=1, budget_gcp=1;
    if (argc > 5) {
        budget = std::stod(argv[2]);
        budget_greedy = std::stod(argv[3]);
        budget_genetic = std::stod(argv[4]);
        budget_gcp = std::stod(argv[5]);
        // is_deepslow = std::stoi(argv[5]);
    }

    std::cout << "Input Parameters:\n";
    std::cout << "  file       = " << file       << "\n";
    std::cout << "  budget     = " << budget     << "\n";
    std::cout << "  slew_rate  = " << slew_rate  << "\n";
    std::cout << "  dwell_time = " << dwell_time << "\n\n";

    test_algorithms (file, out_file, budget, slew_rate, dwell_time, is_deepslow);

    // test_algorithms2(file, out_file, budget, budget_greedy, budget_genetic, budget_gcp, slew_rate, dwell_time, is_deepslow);

    // for(int b=0; b<=60; b+=10){
        // test_algorithms (file, out_file, double(b), slew_rate, dwell_time, is_deepslow);
    // }
    // std::vector<double> budgets = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200};
    // run_multi_deadlines(file, out_file, out_file2, budgets, slew_rate, dwell_time, is_deepslow);
    return 0;
}


