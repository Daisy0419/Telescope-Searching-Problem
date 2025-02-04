#pragma once

#include <vector>
 

 double calculateBound(const std::vector<std::vector<double>>& costs, 
                        const std::vector<double>& prizes,
                        const std::vector<bool>& unvisited, 
                        int current, 
                        double budget);

std::pair<double, std::vector<int>> branchAndBound(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, 
    double budget,
    int current_tile, 
    std::vector<bool>& unvisited, 
    std::vector<int>& cur_path, 
    double cur_prize, 
    std::vector<int>& best_path, 
    double& best_prize);


std::vector<int> run_branch_bound(const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, double budget, int current_tile, std::vector<int> best_path,
    double best_prize = 0);

void branchAndBoundParallelImpl(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes,
    double budget,
    int current_tile,
    std::vector<bool> unvisited, 
    std::vector<int> path, 
    double prize,    
    double &global_best_prize,
    std::vector<int>& global_best_path);


std::pair<double, std::vector<int>> branchAndBoundParallel(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, 
    double budget,
    int current_tile,
    std::vector<bool> unvisited, 
    std::vector<int> path, 
    double prize);

