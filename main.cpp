#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "Christofides.h"
#include "kMST.h"
#include "ReadData.h"

#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

void print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                std::vector<int> best_path) {
    double sum_probability = 0.0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile << " ";
        sum_probability += probability[tile]; // Accumulate probability
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        total_cost += costs[best_path[i]][best_path[i + 1]]; // Compute total cost
    }

    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; // Correct value
    std::cout << "Total Cost: " << total_cost << std::endl;
}


int main() {
    // std::string file = "../data/7dt_combined.csv";
    // std::string file = "../data/7dt_separate.csv";
    std::string file = "../data/deep_slow.csv";
    double budget = 500;
    double slew_rate = 1;
    double dwell_time = 1;

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;

    // read_data(file, costs, probability, slew_rate, dwell_time);
    read_data_deep_slow(file, costs, probability, slew_rate, dwell_time);
    double sum = 0.0;
    for(int i = 0; i < 20; ++i) {
        sum += probability[i];
    }
    std::cout << sum << std::endl;
        
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs);
    ac.ant_colony_optimization(budget, probability);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;



    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    genetic_optimization(costs, probability, budget, 0);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> path = prize_greedy_path(costs, probability, 0, budget);
    end = std::chrono::high_resolution_clock::now();
    print_path(costs, probability, path);
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

