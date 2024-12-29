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


// Function to choose the next city for an ant
std::pair<int, double> AntColony::choose_next_city(int current_city, const std::vector<int> &visited, double budget) {
    std::vector<double> probabilities(NUM_CITIES, 0.0);
    double sum = 0.0;

    for (int j = 0; j < NUM_CITIES; j++) {
        if (distances[current_city][j] <= budget &&
            std::find(visited.begin(), visited.end(), j) == visited.end()) {
            probabilities[j] = std::pow(pheromones[current_city][j], ALPHA) *
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

// Function to simulate a single ant's tour
std::vector<int> AntColony::ant_tour(double budget) {
    // std::default_random_engine generator;
    // std::uniform_int_distribution<int> distribution(0, NUM_CITIES - 1);
    // int start_city = distribution(generator);
    // std::cout << "start_city: " << start_city << std::endl;
    std::default_random_engine generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(0, NUM_CITIES - 1);
    int start_city = distribution(generator);
    // std::cout << "start_city: " << start_city << std::endl;

    std::vector<int> visited;
    visited.push_back(start_city);
    int current_city = start_city;
    double remaining_budget = budget;

    while (remaining_budget > 0) {
        auto [next_city, cost] = choose_next_city(current_city, visited, remaining_budget);
        if (next_city == -1) break;

        visited.push_back(next_city);
        current_city = next_city;
        remaining_budget -= cost;
    }

    return visited;
}


double AntColony::simulate_ants(std::vector<int> &best_tour, double initial_budget) {
    double best_prize = 0;
    std::vector<std::vector<int>> visited(NUM_ANTS);

    // Parallelize the ant tours using OpenMP
    // #pragma omp parallel for
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        visited[ant] = ant_tour(initial_budget);
    }

    // Evaporate pheromones
    // #pragma omp parallel for
    for (int i = 0; i < NUM_CITIES; i++) {
        #pragma omp parallel for
        for (int j = 0; j < NUM_CITIES; j++) {
            pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
        }
    }

    // Update pheromones
    // #pragma omp parallel for reduction(max : best_prize)
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        if (visited[ant].size() < 2) continue;

        double tour_cost = 0.0;
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            tour_cost += distances[visited[ant][i]][visited[ant][i + 1]];
        }

        double tour_prize = visited[ant].size();

        // #pragma omp critical
        if (tour_prize > best_prize) {
            best_prize = tour_prize;
            best_tour = visited[ant];
        }

        double contribution = Q * tour_prize / tour_cost;
        // #pragma omp critical
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
        }
    }

    return best_prize;
}



// double AntColony::simulate_ants(std::vector<int> &best_tour, double initial_budget) {
//     double best_prize = 0;
//     std::vector<std::vector<int>> visited(NUM_ANTS);

//     // Parallelize the ant tours using OpenMP
//     for (int ant = 0; ant < NUM_ANTS; ant++) {
//         visited[ant] = ant_tour(initial_budget);
//     }

//     // Evaporate pheromones
//     for (int i = 0; i < NUM_CITIES; i++) {
//         for (int j = 0; j < NUM_CITIES; j++) {
//             pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
//         }
//     }

//     // Update pheromones
//     for (int ant = 0; ant < NUM_ANTS; ant++) {
//         if (visited[ant].size() < 2) continue;

//         double tour_cost = 0.0;
//         for (size_t i = 0; i < visited[ant].size() - 1; i++) {
//             tour_cost += distances[visited[ant][i]][visited[ant][i + 1]];
//         }

//         double tour_prize = visited[ant].size();


//         if (tour_prize > best_prize) {
//             best_prize = tour_prize;
//             best_tour = visited[ant];
//         }

//         double contribution = Q * tour_prize / tour_cost;

//         for (size_t i = 0; i < visited[ant].size() - 1; i++) {
//             pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
//         }
//     }

//     return best_prize;
// }


void AntColony::ant_colony_optimization(double budget) {
    std::vector<int> best_tour;
    double best_prize = 0;
    // double budget = 47358.8;

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        double prize = simulate_ants(best_tour, budget);
        best_prize = std::max(best_prize, prize);

        // std::cout << "Iteration " << iteration + 1 << " Best Prize: " << best_prize << std::endl;
    //     std::cout << "Best Tour: ";
    //     for (int city : best_tour) {
    //         std::cout << city << " ";
    //     }
    //     std::cout << std::endl;
    }
    std::cout << "Best Prize: " << best_tour.size() << std::endl;
    std::cout << "Best Tour: ";
    for (int city : best_tour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_tour.size() - 1; i++) {
        total_cost += distances[best_tour[i]][best_tour[i + 1]];
    }

    std::cout << "Total Cost: " << total_cost << std::endl;
}


