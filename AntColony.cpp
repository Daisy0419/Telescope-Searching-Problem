#include "AntColony.h"
#include "helplers.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

AntColony::AntColony(std::vector<std::vector<double>>& costs) {
    distances = costs;
    NUM_CITIES = costs.size();
    NUM_ANTS = 2 * NUM_CITIES;

    heuristic.resize(NUM_CITIES, std::vector<double>(NUM_CITIES, 0.0));

    initialize_heuristic();
    pheromones.resize(NUM_CITIES, std::vector<double>(NUM_CITIES, 0.5));
}


void AntColony::initialize_heuristic() {
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            heuristic[i][j] = i != j ? 1.0 / distances[i][j] : 0.0;
        }
    }
}

std::pair<int, double> AntColony::choose_next_city(int current_city, const std::vector<int>& visited, double budget, const std::vector<double>& prizes) {
    std::vector<double> probabilities(NUM_CITIES, 0.0);
    double sum = 0.0;

    for (int j = 0; j < NUM_CITIES; j++) {
        if (distances[current_city][j] <= budget &&
            std::find(visited.begin(), visited.end(), j) == visited.end()) {
            probabilities[j] = std::pow(pheromones[current_city][j], ALPHA) *
                               std::pow(prizes[j], GAMMA) *  // Influence of prize
                               std::pow(1.0 / distances[current_city][j], BETA);
            sum += probabilities[j];
        }
    }

    if (sum == 0.0) return {-1, 0.0}; // No valid moves

    for (int j = 0; j < NUM_CITIES; j++) {
        probabilities[j] /= sum;
    }

    double rand_prob = random_double(0.0, 1.0);
    double cumulative = 0.0;
    for (int j = 0; j < NUM_CITIES; j++) {
        cumulative += probabilities[j];
        if (rand_prob <= cumulative) {
            return {j, distances[current_city][j]};
        }
    }

    return {-1, 0.0};
}

std::vector<int> AntColony::ant_tour(double budget, const std::vector<double>& prizes, int start_city) {
    std::vector<int> visited;
    visited.push_back(start_city);
    int current_city = start_city;
    double remaining_budget = budget;

    while (remaining_budget > 0) {
        auto [next_city, cost] = choose_next_city(current_city, visited, remaining_budget, prizes);
        if (next_city == -1) break;

        visited.push_back(next_city);
        current_city = next_city;
        remaining_budget -= cost;
    }

    return visited;
}

double AntColony::simulate_ants(std::vector<int>& best_tour, double initial_budget, const std::vector<double>& prizes) {
    double best_prize = 0;
    std::vector<std::vector<int>> visited(NUM_ANTS);

    // Simulate each ant's tour
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        visited[ant] = ant_tour(initial_budget, prizes, 0);
    }

    // Evaporate pheromones
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
        }
    }

    // Update pheromones
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        if (visited[ant].size() < 2) continue;

        double tour_cost = 0.0;
        double tour_prize = 0.0;
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            tour_cost += distances[visited[ant][i]][visited[ant][i + 1]];
        }
        for (int city : visited[ant]) {
            tour_prize += prizes[city];
        }

        if (tour_prize > best_prize) {
            best_prize = tour_prize;
            best_tour = visited[ant];
        }

        double contribution = Q * tour_prize / tour_cost;
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
        }
    }

    return best_prize;
}

std::vector<int> AntColony::ant_colony_optimization(double budget, const std::vector<double>& prizes) {
    std::vector<int> best_tour;
    double best_prize = 0;

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        double prize = simulate_ants(best_tour, budget, prizes);
        best_prize = std::max(best_prize, prize);
    }
    return best_tour;


    // std::cout << "Best Prize: " << best_prize << std::endl;
    // std::cout << "Num Tiles in Path: " << best_tour.size() << std::endl;
    // // std::cout << "Num Tiles in Path: " << best_tour.size() << std::endl;
    // std::cout << "Best Tour: ";
    // for (int city : best_tour) {
    //     std::cout << city << " ";
    // }
    // std::cout << std::endl;

    // double total_cost = 0.0;
    // for (size_t i = 0; i < best_tour.size() - 1; i++) {
    //     total_cost += distances[best_tour[i]][best_tour[i + 1]];
    // }

    // std::cout << "Total Cost: " << total_cost << std::endl;
}

