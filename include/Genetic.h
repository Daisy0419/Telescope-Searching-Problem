#pragma once

#include <vector>
#include <string>
#include <unordered_set>

// std::pair<double, double> fitness(const std::vector<int>& path, const std::vector<std::vector<double>>& costs, 
//                                 const std::vector<double>& prizes);


// std::vector<int> genetic_optimization(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
//                           double budget, int start_city, std::vector<std::vector<int>> init_paths={});

// std::vector<int> get_best_path(const std::vector<std::vector<int>>& paths, const std::vector<double>& prizes);

// void my_print_path(const std::vector<std::vector<double>>& distances, const std::vector<double>& probability, 
//                 std::vector<int> best_path); 

// std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, 
//                                 const std::vector<std::vector<double>>& costs, 
//                                 const std::vector<double>& prizes, size_t top_n);

// std::vector<int> get_best_path(const std::vector<std::vector<int>>& paths, const std::vector<double>& prizes);

std::vector<int> genetic_optimization_st(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                          double budget, int start_city, int end_city, std::vector<std::vector<int>> init_paths={});

