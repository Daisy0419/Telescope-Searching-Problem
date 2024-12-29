#pragma once
#include <vector>

std::vector<int> greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget);
std::vector<int> max_unrooted_path(const std::vector<std::vector<double>>& distances, int budget);
void max_rooted_path(const std::vector<std::vector<double>>& distances, int budget, int start_city);
