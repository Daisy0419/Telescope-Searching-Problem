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
        // std::cout << costs[best_path[i]][best_path[i + 1]] << " ";
        total_cost += costs[best_path[i]][best_path[i + 1]]; 
    }
    // std::cout << std::endl;
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
                  std::vector<int> rank,
                  std::vector<int> best_path, 
                  double elapsed_time, int init_pos_idx) {

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
    double total_cost = 0.0;

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


void write_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, double slew_rate, double dwell_time, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> rank,
                  std::vector<int> best_path, 
                  double elapsed_time, std::vector<double>& deadlines) {

    bool file_exists = std::filesystem::exists(filename); 

    std::ofstream outfile(filename, std::ios::app);  

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    //Write header
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,SlewRate,DwellTime,NumTiles,SumProb,TotalCost,TimeSec,Deadlines,Path\n";
    }

    double sum_probability = 0.0;
    double total_cost = 0.0;

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


    for (size_t i = 0; i < deadlines.size(); i++) {
        outfile << deadlines[i];
        if (i < deadlines.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << ",";
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

void test_combined_algorithms (std::string file, std::string out_file, double budget, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<std::vector<double>> region_costs;
    std::vector<double> region_probability;
    std::vector<int> region_ranks;
    double region_sum_probability = 0.99;
    double time_limit = 3600;
    double accu_thr = 0.00001;
    int init_pos_idx = 0;

    // int start_idx = buildGraph(file, costs, probability, region_costs, region_probability, ranks, region_ranks,
    //                             slew_rate, dwell_time, is_deepslow, region_sum_probability, init_pos_idx);

    // int start_idx = buildGraph(file, costs, probability, ranks, slew_rate, dwell_time, is_deepslow, init_pos_idx);
    int start_idx = buildGraph_force_initpos(file, costs, probability, ranks, slew_rate, dwell_time, is_deepslow, init_pos_idx);
    // readOPlibFile(file, costs, probability, ranks);
    std::cout << "region_probability.size(): " << region_probability.size() << std::endl;
    std::cout << "probability.size(): " << probability.size() << std::endl;
    budget -= dwell_time;
    // budget = 1290;
    // int start_idx = 0;

    // probability = region_probability;
    // costs = region_costs;
      
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    
    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prize_ratio_greedy_path(costs, probability, start_idx, budget);

    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, ranks, greedy_path, elapsed_seconds.count(), init_pos_idx);

//   // gurobi with time limit
//     std::cout << "*********gurobi solution with time limit*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> ilp_path = gurobiSolve(costs, probability, start_idx, budget, accu_thr, time_limit);
//     // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget);
//     // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget, genetic_path);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, ilp_path, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "Gurobi", file, budget, slew_rate, dwell_time, costs, probability, ranks, ilp_path, elapsed_seconds.count(), init_pos_idx);


    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int>  genetic_path = genetic_optimization(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, ranks, genetic_path, elapsed_seconds.count(), init_pos_idx);

    // // simulated annealing
    // std::cout << "*********running simulated annealing*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // // genetic_optimization(costs, budget, 0);
    // std::vector<int> sa_path = simulated_annealing_optimization(costs, probability, budget, start_idx);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, ranks, sa_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "SimulatedAnnealing", file, budget, slew_rate, dwell_time, costs, probability, ranks, sa_path, elapsed_seconds.count(), init_pos_idx);

    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs, probability);
    std::vector<int> ac_path = ac.ant_colony_optimization(budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    double best_prize = print_path(costs, probability, ranks, ac_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "AntColony", file, budget, slew_rate, dwell_time, costs, probability, ranks, ac_path, elapsed_seconds.count(), init_pos_idx);

   // mst
    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path = mstNaive(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "mst", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path, elapsed_seconds.count(), init_pos_idx);


    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path3 = mstLemon2(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path3, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "christo_lemon", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path3, elapsed_seconds.count(), init_pos_idx);

  // gurobi with time limit
    std::cout << "*********gurobi solution with time limit*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> ilp_path = gurobiSolve(costs, probability, start_idx, budget, accu_thr, time_limit);
    // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget);
    // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget, genetic_path);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, ilp_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Gurobi", file, budget, slew_rate, dwell_time, costs, probability, ranks, ilp_path, elapsed_seconds.count(), init_pos_idx);

//    // bpc_tsp
//     std::cout << "*********running bpc_tsp*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> bpc_tsp_path = BPC_TSP2(costs, probability, budget-dwell_time, 1);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, bpc_tsp_path, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "BPC-TSP", file, budget, slew_rate, dwell_time, costs, probability, ranks, bpc_tsp_path, elapsed_seconds.count(), init_pos_idx);


//    // kmst
//     std::cout << "*********running kmst*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> kmst_path = findBestTour(costs, probability, budget-dwell_time, 1);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, kmst_path, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "kMST", file, budget, slew_rate, dwell_time, costs, probability, ranks, kmst_path, elapsed_seconds.count(), init_pos_idx);
}


void test_multi_deadline (std::string file, std::string out_file, double budget, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<std::vector<double>> region_costs;
    std::vector<double> region_probability;
    std::vector<int> region_ranks;
    double region_sum_probability = 0.99;
    double time_limit = 3600;
    double accu_thr = 0.00001;
    int init_pos_idx = 0;

    // int start_idx = buildGraph(file, costs, probability, region_costs, region_probability, ranks, region_ranks,
    //                             slew_rate, dwell_time, is_deepslow, region_sum_probability, init_pos_idx);

    int start_idx = buildGraph(file, costs, probability, ranks, slew_rate, dwell_time, is_deepslow, init_pos_idx);
    // int start_idx = buildGraph_force_initpos(file, costs, probability, ranks, slew_rate, dwell_time, is_deepslow, init_pos_idx);
    // readOPlibFile(file, costs, probability, ranks);
    std::cout << "region_probability.size(): " << region_probability.size() << std::endl;
    std::cout << "probability.size(): " << probability.size() << std::endl;
    // budget -= dwell_time;
    // budget = 1290;
    // int start_idx = 0;

    // probability = region_probability;
    // costs = region_costs;
      
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    
    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prize_ratio_greedy_path(costs, probability, start_idx, budget);

    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "Greedy", file, budget, slew_rate, dwell_time, costs, probability, ranks, greedy_path, elapsed_seconds.count(), init_pos_idx);

    // // genetic
    // std::cout << "*********running genetic*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // // genetic_optimization(costs, budget, 0);
    // std::vector<int>  genetic_path = genetic_optimization(costs, probability, budget, start_idx);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, ranks, genetic_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // write_result(out_file, "Genetic", file, budget, slew_rate, dwell_time, costs, probability, ranks, genetic_path, elapsed_seconds.count(), init_pos_idx);

    // mst
    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path0 = mstNaive(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path0, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "mstNaive", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path0, elapsed_seconds.count(), init_pos_idx);


   // mst
    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path = mstLemon(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "mstLemon", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path, elapsed_seconds.count(), init_pos_idx);


    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path3 = mstLemon2(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path3, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "mstLemon2", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path3, elapsed_seconds.count(), init_pos_idx);


    std::cout << "*********running mst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> mst_path4 = mstLemon3(costs, probability, budget, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, mst_path4, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "mstLemon3", file, budget, slew_rate, dwell_time, costs, probability, ranks, mst_path4, elapsed_seconds.count(), init_pos_idx);

//   // gurobi with time limit
//     std::cout << "*********gurobi solution with time limit*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> ilp_path = gurobiSolve(costs, probability, start_idx, budget, accu_thr, time_limit);
//     // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget);
//     // std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget, genetic_path);
//     end = std::chrono::high_resolution_clock::now();
//     print_path(costs, probability, ranks, ilp_path, dwell_time);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     write_result(out_file, "Gurobi", file, budget, slew_rate, dwell_time, costs, probability, ranks, ilp_path, elapsed_seconds.count(), init_pos_idx);

}

void test_multi_deadlines(std::string file, std::string out_file, std::vector<double> budgets, 
                                double slew_rate, double dwell_time, bool is_deepslow=false) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability; 
    std::vector<int> ranks;
    int init_pos_idx = 99;

    if(budgets.size() < 3) {
        std::cerr << "budgets.size() must be >= 3";
        return; 
    }

    double budget1 = budgets[0], budget2 = budgets[1], budget3 = budgets[2];
    int start_idx = buildGraph(file, costs, probability, ranks, slew_rate, dwell_time, is_deepslow, init_pos_idx);
    std::cout << "costs[2][3]: " << costs[2][3];
    std::cout << "costs[17][18]: " << costs[17][18];

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    

    std::cout << "*********running d1 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d1_path = mstLemon2(costs, probability, budget1, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d1_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d1", file, budget1, slew_rate, dwell_time, costs, probability, ranks, d1_path, elapsed_seconds.count(), init_pos_idx);

    std::cout << "*********running d2 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d2_path = mstLemon2(costs, probability, budget2, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d2_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d2", file, budget2, slew_rate, dwell_time, costs, probability, ranks, d2_path, elapsed_seconds.count(), init_pos_idx);

    std::cout << "*********running d3 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> d3_path = mstLemon2(costs, probability, budget3, start_idx);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d3_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d3", file, budget3, slew_rate, dwell_time, costs, probability, ranks, d3_path, elapsed_seconds.count(), init_pos_idx);

    std::vector<std::vector<double>> costs_remaining;
    std::vector<double> prizes_remaining;
    std::vector<int> ranks_mapping;

    std::cout << "*********running d1-d2 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    start_idx = removeNodes(probability, costs, d1_path, costs_remaining, prizes_remaining, ranks_mapping);
    std::vector<int> d1_d2_path_ = mstLemon2(costs_remaining, prizes_remaining, budget2-budget1, start_idx);
    std::vector<int> d1_d2_path = d1_path;
    for (int i = 1; i < d1_d2_path_.size(); ++i) {
        d1_d2_path.push_back(ranks_mapping[d1_d2_path_[i]]);
    }
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d1_d2_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d1_d2", file, budget2-budget1, slew_rate, dwell_time, costs, probability, ranks, d1_d2_path, elapsed_seconds.count(), init_pos_idx);

    std::cout << "*********running d2-d3 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    start_idx = removeNodes(probability, costs, d2_path, costs_remaining, prizes_remaining, ranks_mapping);
    std::vector<int> d2_d3_path_ = mstLemon2(costs_remaining, prizes_remaining, budget3-budget2, start_idx);
    std::vector<int> d2_d3_path = d2_path;
    for (int i = 1; i < d2_d3_path_.size(); ++i) {
        d2_d3_path.push_back(ranks_mapping[d2_d3_path_[i]]);
    }
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d2_d3_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d2_d3", file, budget3-budget2, slew_rate, dwell_time, costs, probability, ranks, d2_d3_path, elapsed_seconds.count(), init_pos_idx);

    std::cout << "*********running d1-d2-d3 path*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    start_idx = removeNodes(probability, costs, d1_d2_path, costs_remaining, prizes_remaining, ranks_mapping);
    std::vector<int> d1d2_d3_path_ = mstLemon2(costs_remaining, prizes_remaining, budget3-budget2, start_idx);
    std::vector<int> d1d2_d3_path = d1_d2_path;
    for (int i = 1; i < d1d2_d3_path_.size(); ++i) {
        d1d2_d3_path.push_back(ranks_mapping[d1d2_d3_path_[i]]);
    }
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ranks, d1d2_d3_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    write_result(out_file, "d1_d2_d3", file, budget3-budget2, slew_rate, dwell_time, costs, probability, ranks, d1d2_d3_path, elapsed_seconds.count(), init_pos_idx);

}


int main (int argc, char** argv) {
    // default paramaters
    std::string out_file = "../results/out.csv";
    // std::string file = "../data_RTSS/filtered_GW191219_163120.fits_slow_deep.csv";
    // std::string file = "../data_RTSS/filtered_GW191219_163120.fits_7dt.csv";
    std::string file = "../data_RTSS/filtered_GW200112_155838.fits_7dt.csv";
    // std::string file = "../data_RTSS/filtered_GW200112_155838.fits_slow_deep.csv";
    double budget = 50;
    double slew_rate = 50;
    double dwell_time = 1;
    bool is_deepslow = false;
    
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

    // test_combined_algorithms (file, out_file, budget, slew_rate, dwell_time, is_deepslow);
    // test_all_algorithms (file, out_file, budget, slew_rate, dwell_time, is_deepslow);

    // for(int b=50; b<=5000; b+=50)
    //     test_multi_deadline (file, out_file, double(b), slew_rate, dwell_time, is_deepslow);


    for(int b=50; b<=5000; b+=50)
        test_multi_deadline (file, out_file, double(b), slew_rate, dwell_time, is_deepslow);
    // std::vector<double> budgets = {500, 1500, 3000};
    // test_multi_deadlines(file, out_file, budgets, slew_rate, dwell_time, is_deepslow);

    return 0;
}


