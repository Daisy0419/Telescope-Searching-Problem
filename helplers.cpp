
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

// void fix_cross(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
//     if (costs.empty() || costs.size() != costs[0].size()) {
//         throw std::invalid_argument("Cost matrix must be non-empty and square.");
//     }
//     if (path.size() > costs.size()) {
//         throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
//     }
//     if (path.size() < 4) {
//         return;  
//     }
//     //2-opt
//     for (size_t i = 0; i < path.size() - 1; ++i) {
//         for (size_t j = i + 2; j < path.size() - 1; ++j) {
//             int a = path[i];     
//             int b = path[i + 1]; 
//             int c = path[j]; 
//             int d = path[j + 1]; 

//             double original_cost = costs[a][b] + costs[c][d];
//             double uncrossed_cost = costs[a][d] + costs[c][b];

//             if (uncrossed_cost < original_cost - 1e-9) {
//                 std::reverse(path.begin() + i + 1, path.begin() + j + 1);
//             }
//         }
//     }
// }