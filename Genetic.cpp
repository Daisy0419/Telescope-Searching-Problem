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

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(unvisited.begin(), unvisited.end(), gen);

    for (int next_city : unvisited) {
        if (costs[current_city][next_city] <= budget) {
            visited.push_back(next_city);
            budget -= costs[current_city][next_city];
            current_city = next_city;
        }
    }

    return visited;
}



std::vector<std::vector<int>> random_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                           double budget, int num_paths){
                                            
    std::vector<std::vector<int>> paths(num_paths);

    // #pragma omp parallel for 
    for (int i = 0; i < num_paths; ++i) {
        paths[i] = generate_random_path(costs, current_city, budget);
    }

    return paths;
}

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


std::vector<int> greedy_path(const std::vector<std::vector<double>>& distances, int current_city, double budget) {
    std::vector<int> visited;
    visited.push_back(current_city);
    int num_cities = distances.size();

    while (budget > 0) {
        double min_dis = std::numeric_limits<double>::infinity();
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (std::find(visited.begin(), visited.end(), i) == visited.end() &&
                distances[current_city][i] < min_dis) {
                next_city = i;
                min_dis = distances[current_city][i];
            }
        }

        if(min_dis > budget) 
            break;

        visited.push_back(next_city);
        budget -= min_dis;
        current_city = next_city;

        // std::cout << "budget: " << budget << std::endl;
    }

    return visited;
}


std::vector<std::vector<int>> partial_greedy_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                            double budget, int num_paths) {

    std::vector<std::vector<int>> paths;
    for(int i = 0; i < num_paths; ++i) {
        std::vector<int> path = std::move(generate_random_path(costs, current_city, budget));
        if(path.size() < 3){
            paths.emplace_back(path);
            continue;
        }
        std::vector<int> greedy_segment = unique_random_ints(1, path.size()-1, 2);
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


void repair_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes,
                 std::vector<int>& path, double budget) {

    std::unordered_set<int> unique_cities;
    std::vector<int> cleaned_path;
    for (int city : path) {
        if (unique_cities.find(city) == unique_cities.end()) {
            unique_cities.insert(city);
            cleaned_path.push_back(city);
        }
    }
    path = std::move(cleaned_path);

    double total_cost = 0.0;

    for (size_t i = 1; i < path.size(); ++i) {
        total_cost += costs[path[i - 1]][path[i]];
    }

    while (total_cost > budget && path.size() > 1) {
        size_t least_prize_index = 1; 
        for (size_t i = 2; i < path.size() - 1; ++i) { 
            if (prizes[path[i]] < prizes[path[least_prize_index]]) {
                least_prize_index = i;
            }
        }

        total_cost -= costs[path[least_prize_index - 1]][path[least_prize_index]];
        if (least_prize_index + 1 < path.size()) {
            total_cost -= costs[path[least_prize_index]][path[least_prize_index + 1]];
            total_cost += costs[path[least_prize_index - 1]][path[least_prize_index + 1]];
        }

        path.erase(path.begin() + least_prize_index);
    }

    int current_city = path.back();
    int num_cities = costs.size();
    std::vector<bool> visited(num_cities, false);
    for (int city : path) {
        visited[city] = true;
    }

    while (true) {
        double max_prize = -1.0;
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (!visited[i] && total_cost + costs[current_city][i] <= budget && prizes[i] > max_prize) {
                next_city = i;
                max_prize = prizes[i];
            }
        }

        if (next_city == -1) {
            break; 
        }

        path.push_back(next_city);
        visited[next_city] = true;
        total_cost += costs[current_city][next_city];
        current_city = next_city;
    }
}


std::vector<int> ordered_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size1 = parent1.size();
    int size2 = parent2.size();
    if (size1 < 3) {
        return size1 < size2 ? parent2 : parent1;  
    }
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



int find_next_edge(int current_city,
                       const std::vector<std::vector<double>>& costs, 
                       const std::vector<double>& prizes,
                       const std::unordered_set<int> unvisited,
                       double budget) {
    int next_city = -1;
    double best_score = -std::numeric_limits<double>::infinity();

    for (auto candidate_city : unvisited) {
        double cost = costs[current_city][candidate_city];
        if (cost <= budget) { 
            double score = prizes[candidate_city] / (cost + 1e-6); //heuristic
            if (score > best_score) {
                next_city = candidate_city;
                best_score = score;
            }
        }
    }

    return next_city;
}

std::vector<int> heuristic_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2,
                                      const std::vector<std::vector<double>>& costs, 
                                      const std::vector<double>& prizes, double budget) {

    std::vector<int> offspring;
    int current_city = parent1[0];
    offspring.push_back(current_city);

    std::unordered_set<int> unvisited(parent1.begin(), parent1.end());
    unvisited.insert(parent2.begin(), parent2.end());
    unvisited.erase(current_city);

    double remaining_budget = budget;

    while (!unvisited.empty() && remaining_budget > 0) {
        int next_city = find_next_edge(current_city, costs, prizes, unvisited, remaining_budget);

        if (next_city == -1) {
            break;
        }

        double cost_to_next_city = costs[current_city][next_city];

        offspring.push_back(next_city);
        remaining_budget -= cost_to_next_city;

        unvisited.erase(next_city);

        current_city = next_city;
    }

    return offspring;
}

std::vector<std::vector<int>> cross_over(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
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
            offspring = std::move(ordered_crossover(parent1, parent2));
        else 
            offspring = std::move( heuristic_crossover(parent1, parent2, costs, prizes, budget));

        repair_path(costs, prizes, offspring, budget);
        // offsprings.push_back(offspring);
        offsprings[i] = std::move(offspring);
    }

    return offsprings;
}


std::vector<int> swap_mutate(std::vector<int> path, int num_swap_pairs) {
    if (path.size() < 3 || num_swap_pairs <= 0) {
        return path;
    }
    for(int j = 0; j < num_swap_pairs; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-1, 2);
        std::iter_swap(path.begin()+swap_elements[0], path.begin()+swap_elements[1]);  
    } 
    return path;
}

std::vector<int> insert_mutate(std::vector<int> path, int num_insertion) {
    if (path.size() < 3 || num_insertion <= 0) {
        return path;
    }
    for(int j = 0; j < num_insertion; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-1, 2);
        int city = path[swap_elements[0]];
        path.erase(path.begin() + swap_elements[0]);
        path.insert(path.begin() + swap_elements[1], city);
    } 

    return path;
}

std::vector<int> reverse_mutate(std::vector<int> path) {
    if (path.size() < 3) {
        return path;
    }
    std::vector<int> range = unique_random_ints(1, path.size()-1, 2);
    int head = std::min(range[0], range[1]);
    int tail = std::max(range[0], range[1]);
    std::reverse(path.begin() + head, path.begin() + tail + 1);

    return path;
}

std::vector<std::vector<int>> mutate(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, double budget, int num_mutation) {
    // std::vector<std::vector<int>> mutate_paths;
    std::vector<std::vector<int>> mutate_paths(num_mutation);
    int num_path = paths.size();

    #pragma omp parallel for schedule(dynamic, 2)
    for(int i = 0; i < num_mutation; ++i) {
        int selected_path = random_int(0, paths.size()-1);
        std::vector<int> path;

        if(i < num_mutation/3)
            path = std::move(swap_mutate(paths[selected_path], 5));
        else if(i >= num_mutation/3 && i < 2*num_mutation/3)
            path = std::move(insert_mutate(paths[selected_path], 5));
        else
            path = std::move(reverse_mutate(paths[selected_path]));

        repair_path(costs, prizes, path, budget);
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

// std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, 
//                                                const std::vector<std::vector<double>>& costs, 
//                                                const std::vector<double>& prizes, size_t top_n) {
//     std::vector<std::pair<std::vector<int>, double>> scored_paths;

//     for (const auto& path : paths) {
//         auto [total_prize, total_cost] = fitness(path, costs, prizes);
//         double score = total_prize;
//         scored_paths.push_back({path, score});
//     }

//     std::sort(scored_paths.begin(), scored_paths.end(),
//               [](const auto& a, const auto& b) {
//                   return a.second > b.second; 
//               });

//     std::vector<std::vector<int>> selected_paths;
//     for (size_t i = 0; i < std::min(top_n, scored_paths.size()); ++i) {
//         selected_paths.push_back(scored_paths[i].first);
//     }

//     return selected_paths;
// }

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
void my_print_path(const std::vector<std::vector<double>>& distances, const std::vector<double>& probability, 
                std::vector<int> best_path) {
    std::cout << " Best Prize: " << best_path.size() << std::endl;
    double sum_probability = 0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile << " ";
        sum_probability += probability[tile];
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        total_cost += distances[best_path[i]][best_path[i + 1]];
        // std::cout << total_cost << " ";
    }
    std::cout << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;
}


std::vector<int> evolution(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                           int start_city, int num_path, double budget, int evolution_itr,
                           std::vector<std::vector<int>>& init_paths) {
    // std::vector<std::vector<int>> paths = std::move(partial_greedy_paths(costs, start_city, budget, num_path));
    int num_random_path = num_path - init_paths.size();
    std::vector<std::vector<int>> paths =  std::move(random_paths(costs, start_city, budget, num_path));
    paths.insert(paths.begin(), init_paths.begin(), init_paths.end());
    int num_mutation = 0.5 * num_path, num_cross = 0.5 * num_path;
    int top_n = num_path * 0.94, num_new_gen = num_path - top_n;

    for (int i = 0; i < evolution_itr; ++i) {
        //Mutate and cross-over paths
        std::vector<std::vector<int>> mutate_paths = std::move(mutate(paths, costs, prizes, budget, num_mutation));
        std::vector<std::vector<int>> cross_over_paths = std::move(cross_over(paths, costs, prizes, budget, num_cross));
        
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
        std::vector<std::vector<int>> new_gen_paths =  std::move(random_paths(costs, start_city, budget, num_new_gen));
        // std::vector<std::vector<int>> new_gen_paths = std::move(partial_greedy_paths(costs, start_city, budget, num_new_gen));
        paths.insert(paths.end(), new_gen_paths.begin(), new_gen_paths.end());

    }

    std::vector<int> best_path = get_best_path(paths, prizes);
    // my_print_path(costs, prizes, best_path);
    return best_path;
}


std::vector<int> genetic_optimization(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                          double budget, int start_city, std::vector<std::vector<int>> init_paths) {
    std::vector<int> best_tour = evolution(costs, prizes, start_city, 100, budget, 1000, init_paths);
    return best_tour;

    // double total_prize = 0.0;
    // double total_cost = 0.0;
    // for (size_t i = 0; i < best_tour.size(); ++i) {
    //     total_prize += prizes[best_tour[i]];
    //     if (i > 0) {
    //         total_cost += costs[best_tour[i - 1]][best_tour[i]];
    //     }
    // }

    // std::cout << "Final Best Prize: " << total_prize << std::endl;
    // std::cout << "Total Cost: " << total_cost << std::endl;
    // std::cout << "Num Tiles in Path: " << best_tour.size() << std::endl;
    // std::cout << "Best Tour: ";
    // for (int city : best_tour) {
    //     std::cout << city << " ";
    // }
    // std::cout << std::endl;

    // my_print_path(costs, prizes, best_tour);
}

