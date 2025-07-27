
#include <iostream>
#include <random>
#include <algorithm>

#include "helplers.h"


double random_double(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

int random_int(int min, int max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

std::vector<int> unique_random_ints(int min, int max, int n) {
    if (n > (max - min + 1)) {
        throw std::invalid_argument("Cannot generate more unique numbers than the range size.");
    }
    std::vector<int> numbers(max - min + 1);
    std::iota(numbers.begin(), numbers.end(), min);

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::shuffle(numbers.begin(), numbers.end(), gen);

    return std::vector<int>(numbers.begin(), numbers.begin() + n);
}


// 2-opt, checks all valid pairs of edges in the path
void fix_cross(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    if (costs.empty() || costs.size() != costs[0].size()) {
        throw std::invalid_argument("Cost matrix must be non-empty and square.");
    }
    if (path.size() > costs.size()) {
        throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
    }
    if (path.size() < 4) { // at least 4 nodes
        return;  
    }

    for (size_t i = 1; i < path.size() - 2; ++i) {
        for (size_t j = i + 1; j < path.size() - 2; ++j) {
            int a = path[i - 1]; 
            int b = path[i];   
            int c = path[j];    
            int d = path[j + 1]; 

            double original_cost = costs[a][b] + costs[c][d];
            double uncrossed_cost = costs[a][c] + costs[b][d];

            if (uncrossed_cost < original_cost - 1e-9) {
                std::reverse(path.begin() + i, path.begin() + j + 1);
            }
        }
    }
}

// 2-opt, ensures that the path remains a valid open path from fixed start to fixed end
void fix_cross_st_path(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    if (costs.empty() || costs.size() != costs[0].size()) {
        throw std::invalid_argument("Cost matrix must be non-empty and square.");
    }
    if (path.size() > costs.size()) {
        throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
    }
    if (path.size() < 4) {
        return;  // nothing to optimize
    }

    // Avoid modifying edges adjacent to start (s) or end (t)
    for (size_t i = 1; i < path.size() - 2; ++i) {
        for (size_t j = i + 1; j < path.size() - 1; ++j) {
            int a = path[i - 1];
            int b = path[i];
            int c = path[j];
            int d = path[j + 1];

            if (i == 1 && j + 1 == path.size() - 1) continue; // affect both s and t
            if (i == 1 && a == path[0]) continue;             // affects s
            if (j + 1 == path.size() - 1 && d == path.back()) continue; // affects t

            double original_cost = costs[a][b] + costs[c][d];
            double uncrossed_cost = costs[a][c] + costs[b][d];

            if (uncrossed_cost < original_cost - 1e-9) {
                std::reverse(path.begin() + i, path.begin() + j + 1);
            }
        }
    }
}


// compute total cost of a given path
double calculate_path_cost(const std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    double total_cost = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        total_cost += costs[path[i - 1]][path[i]];
    }
    return total_cost;
}