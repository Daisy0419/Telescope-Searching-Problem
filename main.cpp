#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "Christofides.h"
#include "kMST.h"
#include "ReadData.h"
#include "BranchBound.h"
#include "ILP_scip.h"
#include "ILP_gurobi.h"

#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                std::vector<int> best_path, double dwell_time) {

    if (best_path.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }
    
    double sum_probability = 0.0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile << " ";
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


int main() {
    // std::string file = "../data/7dt_combined.csv";
    std::string file = "../data/7dt_separate.csv";
    // std::string file = "../data/deep_slow.csv";
    double budget = 30;
    double slew_rate = 50;
    double dwell_time = 1;

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;

    read_data(file, costs, probability, slew_rate, dwell_time);
    // read_data_deep_slow(file, costs, probability, slew_rate, dwell_time);

    // std::vector<std::vector<double>> costs = {
    //     {0, 2, 9, 10, 5, 12, 4, 3},
    //     {2, 0, 7, 6, 8, 5, 3, 4},
    //     {9, 7, 0, 4, 10, 2, 6, 9},
    //     {10, 6, 4, 0, 8, 9, 5, 7},
    //     {5, 8, 10, 8, 0, 7, 2, 4},
    //     {12, 5, 2, 9, 7, 0, 6, 3},
    //     {4, 3, 6, 5, 2, 6, 0, 2},
    //     {3, 4, 9, 7, 4, 3, 2, 0}
    // };

    // std::vector<double> probability = {5.0, 8.0, 10.0, 3.0, 7.0, 6.0, 9.0, 4.0};
    // double budget = 15.0; 
    // int start_tile = 0;
    // double dwell_time = 0;
        
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs);
    std::vector<int> ac_path = ac.ant_colony_optimization(budget-dwell_time, probability);
    end = std::chrono::high_resolution_clock::now();
    double best_prize = print_path(costs, probability, ac_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    std::vector<int> genetic_path = genetic_optimization(costs, probability, budget-dwell_time, 0);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, genetic_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;


    // // branch and bound
    // std::cout << "*********running branch and bound*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // std::vector<int> true_path = run_branch_bound(costs, probability, budget-dwell_time, 0, ac_path, best_prize);
    // end = std::chrono::high_resolution_clock::now();
    // print_path(costs, probability, true_path, dwell_time);
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;


    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> greedy_path = prize_greedy_path(costs, probability, 0, budget-dwell_time);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, greedy_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;


    //gurobi
    std::cout << "*********running gurobi solver*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> ilp_path = gurobiSolve(costs, probability, 0, budget-dwell_time);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, ilp_path, dwell_time);
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    // //kmst
    // std::cout << "*********running kmst*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // findBestTour(costs, budget);
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    return 0;
}


