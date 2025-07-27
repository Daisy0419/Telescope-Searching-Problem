#include "helplers.h"
#include "Greedy.h"
#include "Genetic.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <numeric>
#include <omp.h>

// generate a random path from s-t within budget
std::vector<int> generate_random_st_path(const std::vector<std::vector<double>>& costs, int current_city, 
                                        int end_city, double budget) {
    std::vector<int> visited; 
    visited.push_back(current_city);

    int num_cities = costs.size();

    std::vector<int> unvisited(num_cities);
    std::iota(unvisited.begin(), unvisited.end(), 0);
    unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), current_city), unvisited.end());
    unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), end_city), unvisited.end());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(unvisited.begin(), unvisited.end(), gen);

    for (int next_city : unvisited) {
        if (costs[current_city][next_city] + costs[current_city][end_city] <= budget) {
            visited.push_back(next_city);
            budget -= costs[current_city][next_city];
            current_city = next_city;
        }
    }
    visited.push_back(end_city);
    return visited;
}

// batch generate s-t random paths, each within budget
std::vector<std::vector<int>> random_st_paths(const std::vector<std::vector<double>>& costs, 
                                            int current_city, int end_city,
                                            double budget, int num_paths){
                                            
    std::vector<std::vector<int>> paths(num_paths);

    #pragma omp parallel for schedule(dynamic, 2)
    for (int i = 0; i < num_paths; ++i) {
        paths[i] = generate_random_st_path(costs, current_city, end_city, budget);
    }

    return paths;
}


// generate a path/segment of a path from s-t within budget using greedy-in-cost policy
std::vector<int> greedy_st_path(const std::vector<std::vector<double>>& costs, int start_city, 
                                    int end_city, double budget) {
    std::vector<int> visited;
    visited.push_back(start_city);
    int num_cities = costs.size();
    int current_city = start_city;

    while (budget >= costs[current_city][end_city]) {
        double min_dis = std::numeric_limits<double>::infinity();
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (i != end_city && 
                std::find(visited.begin(), visited.end(), i) == visited.end() &&
                costs[current_city][i] < min_dis &&
                (budget > costs[current_city][i] + costs[i][end_city])) {
                next_city = i;
                min_dis = costs[current_city][i];
            }
        }

        if (next_city == -1) {
            break;
        }

        visited.push_back(next_city);
        budget -= min_dis;
        current_city = next_city;

    }
    visited.push_back(end_city);


    return visited;
}

// batch generate paths
std::vector<std::vector<int>> partial_greedy_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                                    int end_city, double budget, int num_paths) {

    std::vector<std::vector<int>> paths;
    for(int i = 0; i < num_paths; ++i) {
        std::vector<int> path = std::move(generate_random_st_path(costs, current_city, end_city, budget));
        if(path.size() < 4){
            paths.emplace_back(path);
            continue;
        }
        std::vector<int> greedy_segment = unique_random_ints(1, path.size()-2, 2);
        int segment1 = std::min(greedy_segment[0], greedy_segment[1]);
        int segment2 = std::max(greedy_segment[0], greedy_segment[1]);

        std::vector<int> head;

        head = std::vector<int>(path.begin(), path.begin() + segment1+1);
        std::vector<int> tail(path.begin()+ segment2, path.end());

        double cost = calculate_path_cost(head, costs) + calculate_path_cost(tail, costs);
        double remaining_budget = budget - cost;

        std::vector<int> middle = greedy_st_path(costs, path[segment1], path[segment2], remaining_budget);

        head.insert(head.end(), middle.begin()+1, middle.end()-1);
        head.insert(head.end(), tail.begin(), tail.end());

        paths.emplace_back(head);
    }
 
    return paths;
}

//repair a path, ensure it include a node at most once; not exceedding budget; and used the budget to maximum
void repair_st_path(const std::vector<std::vector<double>>& costs,
                 const std::vector<double>& prizes,
                 std::vector<int>& path,
                 double budget) {
    if (path.size() < 2) return; 

    int start_city = path.front();
    int end_city = path.back();

    // remove duplicates while preserving order
    std::unordered_set<int> unique_cities;
    std::vector<int> cleaned_path;
    for (int city : path) {
        if (unique_cities.insert(city).second) {
            cleaned_path.push_back(city);
        }
    }

    //ensure start and end stay fixed
    if (cleaned_path.front() != start_city)
        cleaned_path.insert(cleaned_path.begin(), start_city);
    if (cleaned_path.back() != end_city)
        cleaned_path.push_back(end_city);

    path = std::move(cleaned_path);

    //compute initial total cost
    double total_cost = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        int u = path[i - 1], v = path[i];
        if (u < 0 || u >= (int)costs.size() || v < 0 || v >= (int)costs.size()) {
            std::cerr << "[ERROR] Invalid node index in path: " << u << " or " << v
                    << " (costs size = " << costs.size() << "). Path: ";
            for (int node : path) std::cerr << node << " ";
            std::cerr << std::endl;
            return;  // skip or silently drop invalid path
        }
        total_cost += costs[u][v];
    }


    //prune nodes to fit within budget
    while (total_cost > budget && path.size() > 2) {
        size_t least_prize_index = 1;
        for (size_t i = 1; i < path.size() - 1; ++i) {
            if (prizes[path[i]] < prizes[path[least_prize_index]]) {
                least_prize_index = i;
            }
        }

        // Update cost by removing this node and reconnecting neighbors
        int u = path[least_prize_index - 1];
        int v = path[least_prize_index];
        int w = path[least_prize_index + 1];

        total_cost -= costs[u][v];
        total_cost -= costs[v][w];
        total_cost += costs[u][w];

        path.erase(path.begin() + least_prize_index);
    }

    //try inserting nodes between second last and end to improve reward
    int current_city = path[path.size() - 2];  // node before t
    std::vector<bool> visited(costs.size(), false);
    for (int city : path) visited[city] = true;

    while (true) {
        double max_prize = -1.0;
        int next_city = -1;

        for (int i = 0; i < (int)costs.size(); ++i) {
            if (visited[i]) continue;

            double extra_cost = costs[current_city][i] + costs[i][end_city] - costs[current_city][end_city];
            if (total_cost + extra_cost <= budget && prizes[i] > max_prize) {
                max_prize = prizes[i];
                next_city = i;
            }
        }

        if (next_city == -1) break;

        // Insert next_city before end_city
        path.insert(path.end() - 1, next_city);
        visited[next_city] = true;
        total_cost += costs[current_city][next_city] + costs[next_city][end_city] - costs[current_city][end_city];
        current_city = next_city;
    }
}

std::vector<int> ordered_crossover_st(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size1 = parent1.size();
    std::vector<int> offspring;
    offspring.push_back(parent1[0]); // s

    std::unordered_set<int> used{parent1[0], parent1.back()};
    if(size1 >=4){
        auto seg = unique_random_ints(1, size1 - 2, 2);
        int start = std::min(seg[0], seg[1]), end = std::max(seg[0], seg[1]);

        for (int i = start; i <= end; ++i) {
            offspring.push_back(parent1[i]);
            used.insert(parent1[i]);
        }
    } else {
        for (int i = 1; i < size1-1; ++i) {
            offspring.push_back(parent1[i]);
            used.insert(parent1[i]);
        }
    }

    for (int i = 1; i < (int)parent2.size() - 1; ++i) {
        int city = parent2[i];
        if (!used.count(city)) {
            offspring.push_back(city);
            used.insert(city);
        }
    }

    offspring.push_back(parent1.back()); // t
    return offspring;
}


int find_next_city(int current_city, int end_city,
                   const std::vector<std::vector<double>>& costs, 
                   const std::vector<double>& prizes,
                   const std::unordered_set<int>& unvisited,
                   double budget_remaining_to_insert) {
    int next_city = -1;
    double best_score = -std::numeric_limits<double>::infinity();

    for (int candidate_city : unvisited) {
        double cost_to_candidate = costs[current_city][candidate_city];
        double cost_to_end = costs[candidate_city][end_city];
        double total_cost = cost_to_candidate + cost_to_end;

        if (total_cost <= budget_remaining_to_insert) {
            double score = prizes[candidate_city] * 5 / (cost_to_candidate + 1e-6);  // heuristic
            if (score > best_score) {
                best_score = score;
                next_city = candidate_city;
            }
        }
    }

    return next_city;
}

std::vector<int> heuristic_crossover_st(const std::vector<int>& parent1, const std::vector<int>& parent2,
                                        const std::vector<std::vector<double>>& costs, 
                                        const std::vector<double>& prizes,
                                        double budget) {
    int s = parent1[0];
    int t = parent1.back();

    std::vector<int> offspring{ s };
    std::unordered_set<int> unvisited(parent1.begin() + 1, parent1.end() - 1);
    unvisited.insert(parent2.begin() + 1, parent2.end() - 1);

    double remaining_budget = budget;
    int current_city = s;

    while (!unvisited.empty()) {
        double cost_to_end = costs[current_city][t];

        int next_city = find_next_city(current_city, t, costs, prizes, unvisited,
                                       remaining_budget - cost_to_end);

        if (next_city == -1) {
            break;
        }

        double cost_to_next = costs[current_city][next_city];
        offspring.push_back(next_city);
        remaining_budget -= cost_to_next;
        unvisited.erase(next_city);
        current_city = next_city;
    }

    // Add t if budget allows
    double cost_to_t = costs[current_city][t];
    if (cost_to_t <= remaining_budget) {
        offspring.push_back(t);
    }

    return offspring;
}


std::vector<std::vector<int>> cross_over_st(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                        const std::vector<double>& prizes, double budget, int num_cross) {
    // std::vector<std::vector<int>> offsprings;
    if (paths.size() < 2) {
        return std::vector<std::vector<int>>(1, paths[0]);  
    }
    std::vector<std::vector<int>> offsprings(num_cross);

    #pragma omp parallel for schedule(dynamic, 2)
    for (int i = 0; i < num_cross; ++i) {
        std::vector<int> parents = unique_random_ints(0, paths.size() - 1, 2);
        int parent1_idx = parents[0];
        int parent2_idx = parents[1];

        const auto& parent1 = paths[parent1_idx];
        const auto& parent2 = paths[parent2_idx];

        std::vector<int> offspring;
        if (i < num_cross / 2) 
            offspring = std::move(ordered_crossover_st(parent1, parent2));
        else 
            offspring = std::move( heuristic_crossover_st(parent1, parent2, costs, prizes, budget));

        repair_st_path(costs, prizes, offspring, budget);
        // offsprings.push_back(offspring);
        offsprings[i] = std::move(offspring);
    }

    return offsprings;
}


std::vector<int> swap_mutate_st(std::vector<int> path, int num_swap_pairs) {
    if (path.size() < 4 || num_swap_pairs <= 0) {
        return path;
    }
    for(int j = 0; j < num_swap_pairs; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-2, 2);
        std::iter_swap(path.begin()+swap_elements[0], path.begin()+swap_elements[1]);  
    } 
    return path;
}

std::vector<int> insert_mutate_st(std::vector<int> path, int num_insertion) {
    if (path.size() < 4 || num_insertion <= 0) {
        return path;
    }
    for(int j = 0; j < num_insertion; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-2, 2);
        int city = path[swap_elements[0]];
        path.erase(path.begin() + swap_elements[0]);
        path.insert(path.begin() + swap_elements[1], city);
    } 

    return path;
}

std::vector<int> reverse_mutate_st(std::vector<int> path) {
    if (path.size() < 4) {
        return path;
    }
    std::vector<int> range = unique_random_ints(1, path.size()-2, 2);
    int head = std::min(range[0], range[1]);
    int tail = std::max(range[0], range[1]);
    std::reverse(path.begin() + head, path.begin() + tail + 1);

    return path;
}

std::vector<std::vector<int>> mutate_st(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, double budget, int num_mutation) {
    // std::vector<std::vector<int>> mutate_paths;
    std::vector<std::vector<int>> mutate_paths(num_mutation);
    int num_path = paths.size();

    #pragma omp parallel for schedule(dynamic, 2)
    for(int i = 0; i < num_mutation; ++i) {
        int selected_path = random_int(0, paths.size()-1);
        std::vector<int> path;

        if(i < num_mutation/3)
            path = std::move(swap_mutate_st(paths[selected_path], 5));
        else if(i >= num_mutation/3 && i < 2*num_mutation/3)
            path = std::move(insert_mutate_st(paths[selected_path], 5));
        else
            path = std::move(reverse_mutate_st(paths[selected_path]));

        repair_st_path(costs, prizes, path, budget);
        // mutate_paths.emplace_back(path);
        mutate_paths[i] = std::move(path);
    }
    return mutate_paths;
}


std::pair<double, double> fitness(const std::vector<int>& path, const std::vector<std::vector<double>>& costs, 
                                const std::vector<double>& prizes) {
    double total_prize = 0.0;
    double total_cost = 0.0;

    for (size_t i = 0; i < path.size(); ++i) {
        total_prize += prizes[path[i]];
    }

    for (size_t i = 1; i < path.size(); ++i) {
        total_cost += costs[path[i - 1]][path[i]];
    }

    return {total_prize, total_cost};

}

std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, 
                                               const std::vector<std::vector<double>>& costs, 
                                               const std::vector<double>& prizes, size_t top_n) {
    std::vector<std::tuple<std::vector<int>, double, double>> scored_paths;

    for (const auto& path : paths) {
        auto [total_prize, total_cost] = fitness(path, costs, prizes);
        scored_paths.push_back({path, total_prize, total_cost});
    }

    std::sort(scored_paths.begin(), scored_paths.end(),
              [](const auto& a, const auto& b) {
                    if (std::get<1>(a) == std::get<1>(b)) {
                        return std::get<2>(a) < std::get<2>(b); 
                    }
                  return std::get<1>(a) > std::get<1>(b);
              });

    std::vector<std::vector<int>> selected_paths;
    for (size_t i = 0; i < std::min(top_n, scored_paths.size()); ++i) {
        selected_paths.push_back(std::get<0>(scored_paths[i]));
    }

    return selected_paths;
}

std::vector<int> get_best_path(const std::vector<std::vector<int>>& paths, const std::vector<double>& prizes){

    double best_prize = -1.0;
    std::vector<int> best_path;
    for (const auto& path : paths) {
        double current_prize = 0.0;
        for (int city : path) {
            current_prize += prizes[city];
        }
        if (current_prize > best_prize) {
            best_prize = current_prize;
            best_path = path;
        }
    }
    return best_path;
}

std::vector<int> evolution_st(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                           int start_city, int end_city, int num_path, double budget, int evolution_itr,
                           std::vector<std::vector<int>>& init_paths) {
    // std::vector<std::vector<int>> paths = std::move(partial_greedy_paths(costs, start_city, budget, num_path));
    int num_random_path = num_path - init_paths.size();
    std::vector<std::vector<int>> paths =  std::move(random_st_paths(costs, start_city, end_city, budget, num_random_path));
    paths.insert(paths.begin(), init_paths.begin(), init_paths.end());
    int num_mutation = 0.5 * num_path, num_cross = 0.5 * num_path;
    int top_n = num_path * 0.94, num_new_gen = num_path - top_n;

    for (int i = 0; i < evolution_itr; ++i) {
        //Mutate and cross-over paths
        std::vector<std::vector<int>> mutate_paths = std::move(mutate_st(paths, costs, prizes, budget, num_mutation));
        std::vector<std::vector<int>> cross_over_paths = std::move(cross_over_st(paths, costs, prizes, budget, num_cross));
        
        #pragma omp parallel
        {
            #pragma omp single
            {
                // Process mutate_paths
                for(auto &p : mutate_paths) {
                    #pragma omp task
                    fix_cross(p, costs);
                }
                // Process cross_over_paths
                for(auto &p : cross_over_paths) {
                    #pragma omp task
                    fix_cross(p, costs);
                }
            }
        }

        //Combine all paths
        std::vector<std::vector<int>> all_paths = paths;
        all_paths.insert(all_paths.end(), mutate_paths.begin(), mutate_paths.end());
        all_paths.insert(all_paths.end(), cross_over_paths.begin(), cross_over_paths.end());

        //top paths
        paths = std::move(select_top_paths(all_paths, costs, prizes, top_n));
        std::vector<std::vector<int>> new_gen_paths =  std::move(random_st_paths(costs, start_city, end_city, budget, num_new_gen));
        // std::vector<std::vector<int>> new_gen_paths = std::move(partial_greedy_paths(costs, start_city, budget, num_new_gen));
        paths.insert(paths.end(), new_gen_paths.begin(), new_gen_paths.end());

    }

    std::vector<int> best_path = get_best_path(paths, prizes);
    // my_print_path(costs, prizes, best_path);
    return best_path;
}


std::vector<int> genetic_optimization_st(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                          double budget, int start_city, int end_city, std::vector<std::vector<int>> init_paths) {
    std::vector<int> best_tour = evolution_st(costs, prizes, start_city, end_city, 100, budget, 1000, init_paths);
    return best_tour;

}

