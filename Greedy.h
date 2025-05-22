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

std::vector<int> prizeRatioGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   double budget, int start_tile, int end_tile);

std::vector<int> prizeGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                         const std::vector<double>& probability, 
                                         double budget, int start_tile, int end_tile);