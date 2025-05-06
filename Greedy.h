#pragma once
#include <vector>

std::vector<int> prize_greedy_path(const std::vector<std::vector<double>>& costs,
                                   const std::vector<double>& probability,
                                   int init_tile, double budget);

std::vector<int> partial_prize_greedy_path(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   std::vector<int>& visited,
                                   int init_tile, double budget);

std::vector<int> prize_ratio_greedy_path(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   int init_tile, double budget);

std::vector<int> cost_greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget);
void sortWithIndices(std::vector<double>& values, std::vector<int>& indices);