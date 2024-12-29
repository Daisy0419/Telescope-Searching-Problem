#pragma once

#include <vector>
#include <string>

double calculate_path_cost(const std::vector<int>& path, const std::vector<std::vector<double>>& costs);

std::vector<int> generate_random_path(const std::vector<std::vector<double>>& costs, int current_city, double budget);

std::vector<std::vector<int>> random_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                            double budget, int num_paths);

// std::vector<int> greedy_path(const std::vector<std::vector<double>>& costs, int current_city, double budget);

std::vector<int> greedy_path(const std::vector<std::vector<double>>& costs, int start_city, int end_city, double budget);

std::vector<std::vector<int>> partial_greedy_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                            double budget, int num_paths);

void repair_path(const std::vector<std::vector<double>>& costs, std::vector<int>& path, double budget);

std::vector<int> ordered_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2);

std::vector<int> swap_cross_over(std::vector<int>& parent1, std::vector<int>& parent2);

int find_shortest_edge(const std::vector<int>& current_path, int current_city,
                       const std::vector<std::vector<double>>& costs, const std::vector<int>& parent1, 
                       const std::vector<int>& parent2);

std::vector<int> heuristic_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2,
                                      const std::vector<std::vector<double>>& costs, double budget);

std::vector<std::vector<int>> cross_over(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                        double budget, int num_cross);

std::vector<int> swap_mutate(std::vector<int> path, int num_swap_pairs);

std::vector<int> insert_mutate(std::vector<int> path, int num_insertion);

std::vector<int> reverse_mutate(std::vector<int> path);

std::vector<std::vector<int>> mutate(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    double budget, int num_mutation);

std::pair<int, double> fitness(const std::vector<int>& path, const std::vector<std::vector<double>>& costs);

std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                                size_t top_n);

std::vector<int> evolution(const std::vector<std::vector<double>>& costs, int start_city, int num_path, double budget, int evolution_itr);

void genetic_optimization(const std::vector<std::vector<double>>& costs, double budget, int root);