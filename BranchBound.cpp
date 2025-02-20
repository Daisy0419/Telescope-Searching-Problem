// #include <omp.h>
#include <iostream>
#include <climits>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>

double calculateBound(const std::vector<std::vector<double>>& costs, 
                         const std::vector<double>& prizes,
                         const std::vector<bool>& unvisited, 
                         int current, 
                         double budget) {

    std::vector<int> nodes;
    nodes.push_back(current);
    for (int i = 0; i < unvisited.size(); ++i) {
        if (unvisited[i]) {
            nodes.push_back(i);
        }
    }

    if (nodes.size() <= 1) {
        return prizes[current];
    }

    int N = nodes.size();
    std::vector<bool> inMST(N, false);
    std::vector<double> minEdge(N, std::numeric_limits<double>::infinity());

    minEdge[0] = 0.0; 
    double totalCost = 0.0;
    int countAdded = 0;

    for (int step = 0; step < N; ++step) {
        int u = -1;
        for (int i = 0; i < N; ++i) {
            if (!inMST[i] && (u == -1 || minEdge[i] < minEdge[u])) {
                u = i;
            }
        }
        if (u == -1 || totalCost + minEdge[u] > budget) {
            break;
        }
        inMST[u] = true;
        totalCost += minEdge[u];
        countAdded++;

        for (int v = 0; v < N; ++v) {
            if (!inMST[v]) {
                double c = costs[nodes[u]][nodes[v]];
                if (c < minEdge[v]) {
                    minEdge[v] = c;
                }
            }
        }
    }
    int extraNodes = countAdded - 1; 

    std::vector<double> unvisitedPrizes;
    for (int i = 0; i < (int)unvisited.size(); ++i) {
        if (unvisited[i]) {
            unvisitedPrizes.push_back(prizes[i]);
        }
    }
    std::sort(unvisitedPrizes.begin(), unvisitedPrizes.end(), std::greater<double>());
    double boundPrize = 0.0;
    for (int i = 0; i < extraNodes && i < (int)unvisitedPrizes.size(); ++i) {
        boundPrize += unvisitedPrizes[i];
    }
    boundPrize += prizes[current];
    return boundPrize;
}

std::pair<double, std::vector<int>> branchAndBound(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, 
    double budget,
    int current_tile, 
    std::vector<bool>& unvisited, 
    std::vector<int> cur_path, 
    double cur_prize, 
    std::vector<int>& best_path, 
    double& best_prize) {
    
    std::vector<int> cur_best_path = cur_path;
    double cur_best_prize = cur_prize; 

    for (size_t i = 0; i < unvisited.size(); i++) {
        if (unvisited[i] && costs[current_tile][i] < budget) {
            double new_budget = budget - costs[current_tile][i];
            
            double bound = calculateBound(costs, prizes, unvisited, i, new_budget);
            // std::cout << "best_prize: " << best_prize << std::endl;
            // std::cout << "bound + cur_prize: " << bound + cur_prize << std::endl;
            
            if (bound + cur_prize > best_prize) {
                unvisited[i] = false;
                cur_path.push_back(i);
                cur_prize += prizes[i];
                
                std::pair<double, std::vector<int>> result 
                    = branchAndBound(costs, prizes, new_budget, i,
                                     unvisited,
                                     cur_path,   
                                     cur_prize, 
                                     best_path,
                                     best_prize);
                // std::cout << "prize : " << result.first << std::endl;
                // std::cout << "best prize : " << best_prize << std::endl;
                if (result.first >= best_prize) {
                    best_prize = result.first;
                    best_path = result.second;
                }
                if (result.first > cur_best_prize) {
                    cur_best_prize = result.first;
                    cur_best_path = result.second;
                }

                cur_path.pop_back();
                unvisited[i] = true;
                cur_prize -= prizes[i];
            }
        }
    }
    return {cur_best_prize, cur_best_path};
}


std::vector<int> run_branch_bound(const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, double budget, int current_tile, std::vector<int> best_path,
    double best_prize) {
    std::vector<bool> unvisited(prizes.size(), true);
    unvisited[current_tile] = false;
    
    std::vector<int> cur_path;
    cur_path.push_back(current_tile);
    double cur_prize = prizes[current_tile];


    std::pair<double, std::vector<int>> result = branchAndBound(
        costs, prizes, budget, current_tile,
        unvisited,   
        cur_path,  
        cur_prize,  
        best_path,   
        best_prize 
    );
    std::cout << "final best price: " << result.first << std::endl;
    for(auto t : result.second) {
        std::cout << t << " ";
    }
    std::cout << std::endl;
    return best_path;
}


void branchAndBoundParallelImpl(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes,
    double budget,
    int current_tile,
    std::vector<bool> unvisited, 
    std::vector<int> path, 
    double prize,    
    double &global_best_prize,
    std::vector<int>& global_best_path  
) {
    #pragma omp critical
    {
        if (prize > global_best_prize) {
            global_best_prize = prize;
            global_best_path = path;
        }
    }

    for (size_t i = 0; i < unvisited.size(); i++) {
        if (unvisited[i] && costs[current_tile][i] < budget) {
            double new_budget = budget - costs[current_tile][i];
            double bound = calculateBound(costs, prizes, unvisited, i, new_budget);
            
            if (prize + bound > global_best_prize) {
                auto new_unvisited = unvisited;
                auto new_path = path;
                double new_prize = prize;
                
                new_unvisited[i] = false;
                new_path.push_back(i);
                new_prize += prizes[i];

                #pragma omp task firstprivate(new_unvisited, new_path, new_prize, new_budget, i) \
                                 shared(global_best_prize, global_best_path, costs, prizes)
                {
                    branchAndBoundParallelImpl(costs, prizes, new_budget, i,
                                               new_unvisited, new_path, new_prize,
                                               global_best_prize, global_best_path);
                }
            }
        }
    }
    #pragma omp taskwait
}


std::pair<double, std::vector<int>> branchAndBoundParallel(
    const std::vector<std::vector<double>>& costs, 
    const std::vector<double>& prizes, 
    double budget,
    int current_tile,
    std::vector<bool> unvisited, 
    std::vector<int> path, 
    double prize) {

    double global_best_prize = prize;
    std::vector<int> global_best_path = path;

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            branchAndBoundParallelImpl(costs, prizes, budget, current_tile,
                                       unvisited, path, prize,
                                       global_best_prize, global_best_path);
        }
    }
    return {global_best_prize, global_best_path};
}

// std::vector<int> run_branch_bound_parallel(const std::vector<std::vector<double>>& costs, 
//                                            const std::vector<double>& prizes, 
//                                            double budget, 
//                                            int current_tile) {

//     std::vector<bool> unvisited(costs.size(), true);
//     unvisited[current_tile] = false;

//     std::vector<int> initial_path;
//     double initial_prize = 0.0;

//     std::pair<double, std::vector<int>> result = branchAndBoundParallel(
//         costs, prizes, budget, current_tile, unvisited, initial_path, initial_prize
//     );

//     return result.second;
// }

