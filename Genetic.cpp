#include "helplers.h"
#include "Greedy.h"

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




double calculate_path_cost(const std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    double total_cost = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        total_cost += costs[path[i - 1]][path[i]];
    }
    return total_cost;
}


std::vector<int> generate_random_path(const std::vector<std::vector<double>>& costs, int current_city, double budget) {
    std::vector<int> visited;
    visited.push_back(current_city);
    int num_cities = costs.size();

    std::vector<int> unvisited(num_cities);
    std::iota(unvisited.begin(), unvisited.end(), 0);
    unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), current_city), unvisited.end());

    while (!unvisited.empty() && budget > 0) {
        // Shuffle unvisited cities to randomize selection
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(unvisited.begin(), unvisited.end(), gen);

        bool found_next_city = false;
        for (auto it = unvisited.begin(); it != unvisited.end(); ++it) {
            int next_city = *it;
            // if (current_city >= costs.size() || next_city >= costs[current_city].size()) {
            //     throw std::runtime_error("Index out of bounds");
            // }
            if (costs[current_city][next_city] <= budget) {
                visited.push_back(next_city);
                budget -= costs[current_city][next_city];
                current_city = next_city;

                unvisited.erase(it);
                found_next_city = true;
                break;
            }
        }

        if (!found_next_city) {
            break;
        }
    }

    return visited;
}


std::vector<std::vector<int>> random_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                            double budget, int num_paths) {

    std::vector<std::vector<int>> paths;
    for(int i = 0; i < num_paths; ++i) {
        std::vector<int> path = std::move(generate_random_path(costs, current_city, budget));
        paths.emplace_back(path);
    }
    return paths;
}

// std::vector<int> greedy_path(const std::vector<std::vector<double>>& costs, int current_city, double budget) {
//     std::vector<int> visited;
//     visited.push_back(current_city);
//     int num_cities = costs.size();

//     while (budget > 0) {
//         double min_dis = std::numeric_limits<double>::infinity();
//         int next_city = -1;

//         for (int i = 0; i < num_cities; ++i) {
//             if (std::find(visited.begin(), visited.end(), i) == visited.end() &&
//                 costs[current_city][i] < min_dis) {
//                 next_city = i;
//                 min_dis = costs[current_city][i];
//             }
//         }

//         if(min_dis > budget) 
//             break;

//         visited.push_back(next_city);
//         budget -= min_dis;
//         current_city = next_city;

//         // std::cout << "budget: " << budget << std::endl;
//     }

//     return visited;
// }


std::vector<int> greedy_path(const std::vector<std::vector<double>>& costs, int start_city, int end_city, double budget) {
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

std::vector<std::vector<int>> partial_greedy_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                            double budget, int num_paths) {

    std::vector<std::vector<int>> paths;
    for(int i = 0; i < num_paths; ++i) {
        std::vector<int> path = std::move(generate_random_path(costs, current_city, budget));
        std::vector<int> greedy_segment = unique_random_ints(0, path.size(), 2);
        int segment1 = std::min(greedy_segment[0], greedy_segment[1]);
        int segment2 = std::max(greedy_segment[0], greedy_segment[1]);

        std::vector<int> head;
        if(segment2 != path.size()) {
            head = std::vector<int>(path.begin(), path.begin() + segment1+1);
            std::vector<int> tail(path.begin()+ segment2, path.end());

            double cost = calculate_path_cost(head, costs) + calculate_path_cost(tail, costs);
            double remaining_budget = budget - cost;

            std::vector<int> middle = greedy_path(costs, path[segment1], path[segment2], remaining_budget);

            head.insert(head.end(), middle.begin()+1, middle.end()-1);
            head.insert(head.end(), tail.begin(), tail.end());
        }
        else {
            head = std::vector<int>(path.begin(), path.begin() + segment1+1);
            double cost = calculate_path_cost(head, costs);
            double remaining_budget = budget - cost;

            std::vector<int> middle = greedy_path(costs, path[segment1], remaining_budget);
                        head.insert(head.end(), middle.begin()+1, middle.end());
        }


        paths.emplace_back(head);
    }
 
    return paths;
}


void repair_path(const std::vector<std::vector<double>>& costs, std::vector<int>& path, double budget) {
    double total_cost = 0.0;
    bool to_extend = true;

    // Truncate the path if it exceeds the budget
    for (size_t i = 1; i < path.size(); ++i) {
        // std::cout << "city: " << path[i] << " ";
        if (path[i - 1] >= costs.size() || path[i] >= costs[path[i - 1]].size()) {
            std::cerr << "Invalid index in path: path[i-1] = " << path[i - 1]
                    << ", path[i] = " << path[i] << std::endl;
            throw std::runtime_error("Invalid index in path");
        }
        total_cost += costs[path[i - 1]][path[i]];
        if (total_cost > budget) {
            path.resize(i); 
            to_extend = false;
            break;
        }
    }

    //extend the path using a greedy approach if budget allows
    if(to_extend) {
        int current_city = path.back();
        int num_cities = costs.size();
        std::vector<bool> visited(num_cities, false);
        for (int city : path) {
            visited[city] = true; 
        }

        while (true) {
            double min_cost = std::numeric_limits<double>::infinity();
            int next_city = -1;

            for (int i = 0; i < num_cities; ++i) {
                if (!visited[i] && costs[current_city][i] < min_cost && total_cost + costs[current_city][i] <= budget) {
                    next_city = i;
                    min_cost = costs[current_city][i];
                }
            }

            if (next_city == -1) {
                break;
            }

            path.push_back(next_city);
            visited[next_city] = true;
            total_cost += min_cost;
            current_city = next_city;
        }

    }

}

std::vector<int> ordered_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size1 = parent1.size();
    int size2 = parent2.size();
    std::vector<int> offspring(std::max(size1, size2), -1);

    std::vector<int> selected_segment = unique_random_ints(0, size1 - 1, 2);
    int start = selected_segment[0], end = selected_segment[1];
    if (start > end) std::swap(start, end);

    for (int i = start; i <= end; ++i) {
        offspring[i] = parent1[i];
    }

    int idx = 0;
    for (int i = 0; i < size2; ++i) {
        if (idx == start) {
            idx = end + 1;
        }

        int city = parent2[i];
        if (std::find(offspring.begin(), offspring.end(), city) == offspring.end()) {
            if (idx >= offspring.size()) {
                offspring.resize(idx + 1, -1);
            }
            offspring[idx] = city;
            idx += 1;
        }
    }

    auto it = std::find(offspring.begin(), offspring.end(), -1);
    offspring.resize(std::distance(offspring.begin(), it));
    

    return offspring;
}


std::vector<int> swap_cross_over(std::vector<int>& parent1, std::vector<int>& parent2) {
    //to-do

    return std::vector<int>();
}

int find_shortest_edge(const std::vector<int>& current_path, int current_city,
                       const std::vector<std::vector<double>>& costs, const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int next_city = -1;
    double min_cost = std::numeric_limits<double>::infinity();

    for (const auto& parent : {parent1, parent2}) {
        auto it = std::find(parent.begin(), parent.end(), current_city);
        if (it != parent.end()) {
            int parent_idx = it - parent.begin();

            if (parent_idx + 1 < parent.size()) {
                int candidate_city = parent[parent_idx + 1];
                if (current_city >= costs.size() || candidate_city >= costs[current_city].size()) {
                    std::cerr << "Error: Invalid index in costs. current_city = " << current_city
                            << ", candidate_city = " << candidate_city << std::endl;
                    continue;
                }

                if (std::find(current_path.begin(), current_path.end(), candidate_city) == current_path.end() &&
                    costs[current_city][candidate_city] < min_cost) {
                    next_city = candidate_city;
                    min_cost = costs[current_city][candidate_city];
                }
            }
        }
    }

    return next_city;
}


std::vector<int> heuristic_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2,
                                      const std::vector<std::vector<double>>& costs, double budget) {
    int num_cities = costs.size();
    std::vector<int> offspring;

    std::vector<int> unvisited(num_cities);
    std::iota(unvisited.begin(), unvisited.end(), 0);

    int current_city = parent1[0];
    offspring.push_back(current_city);

    // unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), current_city), unvisited.end());
    auto it = std::remove(unvisited.begin(), unvisited.end(), current_city);
    if (it != unvisited.end()) {
        unvisited.erase(it, unvisited.end());
    }


    while (offspring.size() < parent1.size() + parent2.size()) {
        int next_city = find_shortest_edge(offspring, current_city, costs, parent1, parent2);

        if (next_city == -1) {
            int random_idx = random_int(0, unvisited.size() - 1);
            next_city = unvisited[random_idx];
        }

        offspring.push_back(next_city);

        unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), next_city), unvisited.end());

        current_city = next_city;
    }

    return offspring;
}


std::vector<std::vector<int>> cross_over(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                        double budget, int num_cross) {
    std::vector<std::vector<int>> offsprings;

    for (int i = 0; i < num_cross; ++i) {
        std::vector<int> parents = unique_random_ints(0, paths.size() - 1, 2);
        int parent1_idx = parents[0];
        int parent2_idx = parents[1];

        const auto& parent1 = paths[parent1_idx];
        const auto& parent2 = paths[parent2_idx];

        std::vector<int> offspring;
        if (i < num_cross / 2) 
            offspring = std::move(ordered_crossover(parent1, parent2));
        else 
            offspring = std::move(heuristic_crossover(parent1, parent2, costs, budget));

        repair_path(costs, offspring, budget);
        offsprings.push_back(offspring);

    }

    return offsprings;
}


std::vector<int> swap_mutate(std::vector<int> path, int num_swap_pairs) {
    if (path.size() < 2)
        return path;
    for(int j = 0; j < num_swap_pairs; ++j) {
        std::vector<int> swap_elements = unique_random_ints(0, path.size()-1, 2);
        std::iter_swap(path.begin()+swap_elements[0], path.begin()+swap_elements[1]);  
    } 

    return path;
}

std::vector<int> insert_mutate(std::vector<int> path, int num_insertion) {

    for(int j = 0; j < num_insertion; ++j) {
        std::vector<int> swap_elements = unique_random_ints(0, path.size()-1, 2);
        int city = path[swap_elements[0]];
        path.erase(path.begin() + swap_elements[0]);
        path.insert(path.begin() + swap_elements[1], city);
    } 

    return path;
}

std::vector<int> reverse_mutate(std::vector<int> path) {

    std::vector<int> range = unique_random_ints(0, path.size()-1, 2);
    int head = std::min(range[0], range[1]);
    int tail = std::max(range[0], range[1]);
    std::reverse(path.begin() + head, path.begin() + tail + 1);

    return path;
}

std::vector<std::vector<int>> mutate(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    double budget, int num_mutation) {
    std::vector<std::vector<int>> mutate_paths;
    int num_path = paths.size();

    for(int i = 0; i < num_mutation; ++i) {
        int selected_path = random_int(0, paths.size()-1);
        std::vector<int> path;

        if(i < num_mutation/3)
            path = std::move(swap_mutate(paths[selected_path], 5));
        else if(i >= num_mutation/3 && i < 2*num_mutation/3)
            path = std::move(insert_mutate(paths[selected_path], 5));
        else
            path = std::move(reverse_mutate(paths[selected_path]));

        repair_path(costs, path, budget);
        mutate_paths.emplace_back(path);
    }
    return mutate_paths;
}

std::pair<int, double> fitness(const std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    int visited_cities = path.size();               
    double total_cost = calculate_path_cost(path, costs); 
    return {visited_cities, -total_cost};           
}


std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                                size_t top_n) {
    std::vector<std::pair<std::vector<int>, std::pair<int, double>>> ranked_paths;
    for (const auto& path : paths) {
        ranked_paths.push_back({path, fitness(path, costs)});
    }

    if (ranked_paths.size() > top_n) {
        std::partial_sort(ranked_paths.begin(), ranked_paths.begin() + top_n, ranked_paths.end(),
                          [](const auto& a, const auto& b) {
                              if (a.second.first == b.second.first) {
                                  return a.second.second < b.second.second; 
                              }
                              return a.second.first > b.second.first; 
                          });
    } else {
        std::sort(ranked_paths.begin(), ranked_paths.end(),
                  [](const auto& a, const auto& b) {
                      if (a.second.first == b.second.first) {
                          return a.second.second < b.second.second; 
                      }
                      return a.second.first > b.second.first; 
                  });
    }

    std::vector<std::vector<int>> top_paths;
    for (size_t i = 0; i < std::min(top_n, ranked_paths.size()); ++i) {
        top_paths.push_back(ranked_paths[i].first);
    }

    return top_paths;
}


std::vector<int> evolution(const std::vector<std::vector<double>>& costs, int start_city, int num_path, double budget, int evolution_itr) {
    // std::vector<std::vector<int>> paths = std::move(random_paths(costs, start_city, budget, num_path));
    std::vector<std::vector<int>> paths = std::move(partial_greedy_paths(costs, start_city, budget, num_path));

    // for(auto path : paths) {
    //     std::cout << "path length: " << path.size() << std::endl;
    // }
    paths[num_path-1] = greedy_path(costs, start_city, budget);
    // std::cout << "init path length: " << paths[num_path-1].size() << std::endl;
    std::sort(paths.begin(), paths.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a.size() > b.size();
    });

    int init_max_prize = paths[0].size();

    int num_mutation = 0.4*num_path, num_cross = 0.4*num_path, num_new_gen = 0.2*num_path;
    int top_n = num_path * 0.94, new_random = num_path * 0.02, new_partial_greedy = num_path * 0.02;

    for(int i = 0; i < evolution_itr; ++i) {
        // std::cout << "itr: " << i << std::endl;
        std::vector<std::vector<int>> mutate_paths = std::move(mutate(paths, costs, budget, num_mutation));
        std::vector<std::vector<int>> cross_over_paths = std::move(cross_over(paths, costs, budget, num_cross));
        std::vector<std::vector<int>> new_gen_paths(num_new_gen, std::vector<int>());
        for(int i = 0; i < num_new_gen; ++i) {
            new_gen_paths[i] = generate_random_path(costs, start_city, budget);
        }

        std::vector<std::vector<int>> all_paths;
        all_paths.insert(all_paths.end(), paths.begin(), paths.end());
        all_paths.insert(all_paths.end(), mutate_paths.begin(), mutate_paths.end());
        all_paths.insert(all_paths.end(), cross_over_paths.begin(), cross_over_paths.end());
        all_paths.insert(all_paths.end(), new_gen_paths.begin(), new_gen_paths.end());
       
        paths = std::move(select_top_paths(all_paths, costs, top_n));
        std::vector<std::vector<int>> new_random_paths = std::move(random_paths(costs, start_city, budget, new_random));
        std::vector<std::vector<int>> new_greedy_paths = std::move(partial_greedy_paths(costs, start_city, budget, new_partial_greedy));
        paths.insert(paths.end(), new_random_paths.begin(), new_random_paths.end());
        // paths.insert(paths.end(), new_greedy_paths.begin(), new_greedy_paths.end());
        
        // std::cout << "best Prize: " << paths[0].size() << std::endl;
    }
    std::cout << " Init Best Prize: " << init_max_prize << std::endl;
    return paths[0];
}

void genetic_optimization(const std::vector<std::vector<double>>& costs, double budget, int root){
    std::vector<int> best_tour = evolution(costs, 0, 100, budget, 1000);

    std::cout << " Final Best Prize: " << best_tour.size() << std::endl;
    std::cout << "Best Tour: ";
    for (int city : best_tour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_tour.size() - 1; i++) {
        total_cost += costs[best_tour[i]][best_tour[i + 1]];
        // std::cout << total_cost << " ";
    }
    std::cout << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;
}