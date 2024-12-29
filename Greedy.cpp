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

std::vector<int> greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget) {
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



std::vector<int> max_unrooted_path(const std::vector<std::vector<double>>& distances, int budget) {
    int num_cities = distances.size();
    int prize = 0;
    std::vector<int> path;
    for(int i = 0; i < num_cities; ++i) {
        std::vector<int> cur_path = greedy_path(distances, i, budget);
        if(cur_path.size() > prize){
            path = cur_path;
            prize = cur_path.size();
        }
    }
    return path;
}


void max_rooted_path(const std::vector<std::vector<double>>& distances, int budget, int start_city) {

    std::vector<int> best_tour = greedy_path(distances, start_city, budget);

    std::cout << " Best Prize: " << best_tour.size() << std::endl;
    std::cout << "Best Tour: ";
    for (int city : best_tour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_tour.size() - 1; i++) {
        total_cost += distances[best_tour[i]][best_tour[i + 1]];
        // std::cout << total_cost << " ";
    }
    std::cout << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;
}
