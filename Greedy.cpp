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

void sortWithIndices(std::vector<double>& values, std::vector<int>& indices) {
    std::vector<std::pair<double, int>> valueIndexPairs;
    for (size_t i = 0; i < values.size(); ++i) {
        valueIndexPairs.emplace_back(values[i], i);
    }

    std::sort(valueIndexPairs.begin(), valueIndexPairs.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first;
              });

    for (size_t i = 0; i < valueIndexPairs.size(); ++i) {
        values[i] = valueIndexPairs[i].first;
        indices[i] = valueIndexPairs[i].second;
    }
}

std::vector<int> prize_greedy_path(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   int init_tile, double budget) {

    // Create a sorted list of indices based on prize-to-cost ratio
    std::vector<int> indices(probability.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<std::pair<double, int>> prizeToCost(probability.size());

    for (int i = 0; i < probability.size(); ++i) {
        if (i != init_tile) { // Avoid dividing by zero for self-loop
            double min_cost = *std::min_element(costs[init_tile].begin(), costs[init_tile].end());
            prizeToCost[i] = {probability[i] / std::max(1e-9, min_cost), i};
        } else {
            prizeToCost[i] = {0.0, i}; // Initialize the starting node with a low ratio
        }
    }

    std::sort(prizeToCost.begin(), prizeToCost.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first; // Sort by descending prize-to-cost ratio
              });

    // Greedy path selection
    std::vector<int> path;
    if (init_tile == -1) {
        init_tile = prizeToCost.front().second; // Start with the best initial tile
    }

    int current_tile = init_tile;
    path.push_back(current_tile);

    for (const auto& [_, tile] : prizeToCost) {
        if (tile == current_tile || costs[current_tile][tile] > budget) {
            continue;
        }

        // Add the tile to the path
        budget -= costs[current_tile][tile];
        current_tile = tile;
        path.push_back(current_tile);

        // Stop if no more budget remains
        if (budget <= 0) break;
    }

    return path;
}




std::vector<int> cost_greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget) {
    std::vector<int> visited;
    visited.push_back(current_city);
    int num_cities = distances.size();

    while (budget > 0) {
        double min_dis = std::numeric_limits<double>::infinity();
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (std::find(visited.begin(), visited.end(), i) == visited.end() &&
                distances[current_city][i] < min_dis) {
                next_city = i;
                min_dis = distances[current_city][i];
            }
        }

        if(min_dis > budget) 
            break;

        visited.push_back(next_city);
        budget -= min_dis;
        current_city = next_city;

        // std::cout << "budget: " << budget << std::endl;
    }

    return visited;
}



// std::vector<int> max_unrooted_path(const std::vector<std::vector<double>>& distances, int budget) {
//     int num_cities = distances.size();
//     int prize = 0;
//     std::vector<int> path;
//     for(int i = 0; i < num_cities; ++i) {
//         std::vector<int> cur_path = greedy_path(distances, i, budget);
//         if(cur_path.size() > prize){
//             path = cur_path;
//             prize = cur_path.size();
//         }
//     }
//     return path;
// }