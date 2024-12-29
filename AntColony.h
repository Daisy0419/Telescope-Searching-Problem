#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <utility>

class AntColony {
private:
    int NUM_CITIES = 100;
    int NUM_ANTS = 200;
    int MAX_ITERATIONS = 200;
    double ALPHA = 1.0; // Pheromone influence
    double BETA = 1.0;  // Heuristic influence
    double EVAPORATION_RATE = 0.1;
    double Q = 1.0;   // Pheromone deposit factor
    double V = 10.0; // Telescope Velocity

    std::vector<std::vector<double>> distances;
    std::vector<std::vector<double>> pheromones;
    std::vector<std::vector<double>> heuristic;

    void initialize_heuristic();

    std::pair<int, double> choose_next_city(int current_city, const std::vector<int> &visited, double budget);
    std::vector<int> ant_tour(double budget);
    double simulate_ants(std::vector<int> &best_tour, double initial_budget);

public:
    AntColony(std::vector<std::vector<double>>&);

    void ant_colony_optimization(double);
};
