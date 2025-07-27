#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <utility>

class AntColony {
private:
    int NUM_CITIES;
    int NUM_ANTS = 100;
    int MAX_ITERATIONS = 1000;
    double ALPHA = 1.0; // Pheromone influence
    double BETA = 1.0;  // Heuristic influence
    double GAMMA = 1.0;  // Prize influence
    double EVAPORATION_RATE = 0.1;
    double Q = 20.0;   // Pheromone deposit factor

    std::vector<std::vector<double>> distances;
    std::vector<std::vector<double>> pheromones;
    std::vector<double> prizes;
    std::vector<std::vector<double>> heuristic;

    void initialize_heuristic();

    // for s-t path
    std::pair<int, double> choose_next_city(int current_city, int end_city, 
                                            const std::vector<bool>& visited_set, double budget);
    std::vector<int> ant_tour(double budget, int start_city, int end_city);
    double simulate_ants(std::vector<int> &best_tour, double &best_prize, 
                        double initial_budget, int start_city, int end_city);

public:
    AntColony(std::vector<std::vector<double>> cost, std::vector<double> prizes);
    std::vector<int> ant_colony_optimization(double budget, int start_city, int end_city);
};
