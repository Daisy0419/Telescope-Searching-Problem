#pragma once
#include <vector>

// std::vector<int> greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget);
// std::vector<int> max_unrooted_path(const std::vector<std::vector<double>>& distances, int budget);
// void max_rooted_path(const std::vector<std::vector<double>>& distances, int budget, int start_city);

std::vector<int> prize_greedy_path(const std::vector<std::vector<double>>& costs,
                                   const std::vector<double>& probability,
                                   int init_tile, double budget);

std::vector<int> cost_greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget);
void sortWithIndices(std::vector<double>& values, std::vector<int>& indices);