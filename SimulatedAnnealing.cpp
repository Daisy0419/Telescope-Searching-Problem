#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>

#include "SimulatedAnnealing.h"

std::mt19937 rng(std::random_device{}());

double path_cost(const std::vector<std::vector<double>>& costs, const std::vector<int>& path) {
    double cost_sum = 0.0;
    for (size_t i = 0; i < path.size() - 1; i++)
        cost_sum += costs[path[i]][path[i + 1]];
    return cost_sum;
}

double path_prize(const std::vector<double>& prizes, const std::vector<int>& path) {
    double prize_sum = 0.0;
    for (int city : path)
        prize_sum += prizes[city];
    return prize_sum;
}

std::vector<int> initial_solution(const std::vector<std::vector<double>> &costs, 
                                    int budget, int start_city) {
    int NUM_CITIES = costs.size();               
    std::vector<int> cities(NUM_CITIES);
    std::iota(cities.begin(), cities.end(), 0);
    cities.erase(cities.begin() + start_city); // Remove start city
    std::shuffle(cities.begin(), cities.end(), rng);
    
    std::vector<int> solution = {start_city};
    double current_budget = 0.0;

    for (int city : cities) {
        double extra_cost = costs[solution.back()][city];
        if (current_budget + extra_cost > budget)
            break;
        solution.push_back(city);
        current_budget += extra_cost;
    }
    return solution;
}

std::vector<int> two_opt(const std::vector<int>& path) {
    if (path.size() <= 3) return path;  // No meaningful swap

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 
    
    size_t a = dist(rng);
    size_t b = dist(rng);
    if (a > b) std::swap(a, b);

    std::reverse(new_path.begin() + a, new_path.begin() + b + 1);
    return new_path;
}

std::vector<int> three_opt(const std::vector<int>& path) {
    size_t n = path.size();
    if (n <= 5) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, n - 2);  

    size_t a, b, c;
    do {
        a = dist(rng);
        b = dist(rng);
        c = dist(rng);
    } while (a >= b || b >= c);

    std::vector<int> segment1(new_path.begin(), new_path.begin() + a);
    std::vector<int> segment2(new_path.begin() + a, new_path.begin() + b);
    std::vector<int> segment3(new_path.begin() + b, new_path.begin() + c);
    std::vector<int> segment4(new_path.begin() + c, new_path.end());

    std::uniform_int_distribution<int> case_dist(1, 4);
    int case_num = case_dist(rng);

    switch (case_num) {
        case 1: std::reverse(segment2.begin(), segment2.end()); break;
        case 2: std::reverse(segment3.begin(), segment3.end()); break;
        case 3: std::swap(segment2, segment3); break;
        case 4:
            std::reverse(segment2.begin(), segment2.end());
            std::reverse(segment3.begin(), segment3.end());
            break;
    }

    std::vector<int> candidate;
    candidate.reserve(n);
    candidate.insert(candidate.end(), segment1.begin(), segment1.end());
    candidate.insert(candidate.end(), segment2.begin(), segment2.end());
    candidate.insert(candidate.end(), segment3.begin(), segment3.end());
    candidate.insert(candidate.end(), segment4.begin(), segment4.end());

    return candidate;
}

std::vector<int> swap_cities(const std::vector<int>& path) {
    if (path.size() <= 2) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 1); 

    size_t a = dist(rng);
    size_t b = dist(rng);
    std::swap(new_path[a], new_path[b]);

    return new_path;
}

std::vector<int> shift_segment(const std::vector<int>& path) {
    if (path.size() <= 3) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 

    size_t start = dist(rng);
    size_t end = dist(rng);
    if (start > end) std::swap(start, end);

    std::vector<int> segment(new_path.begin() + start, new_path.begin() + end + 1);
    new_path.erase(new_path.begin() + start, new_path.begin() + end + 1);

    std::uniform_int_distribution<size_t> insert_dist(1, new_path.size()); 
    size_t insert_pos = insert_dist(rng);
    new_path.insert(new_path.begin() + insert_pos, segment.begin(), segment.end());

    return new_path;
}

std::vector<int> shuffle_segment(const std::vector<int>& path) {
    if (path.size() <= 3) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 

    size_t start = dist(rng);
    size_t end = dist(rng);
    if (start > end) std::swap(start, end);

    std::shuffle(new_path.begin() + start, new_path.begin() + end + 1, rng);

    return new_path;
}

std::vector<int> change_city(const std::vector<int>& path,
                             const std::vector<std::vector<double>>& costs,
                             const std::vector<double>& prizes,
                             double budget,
                             std::vector<bool>& visited) {
    if (path.size() < 2) 
        return path;

    std::vector<int> new_path = path;
    double current_cost = path_cost(costs, path);

    // Remove a random city (not the start city at index 0)
    std::uniform_int_distribution<size_t> dist(1, path.size() - 1);
    size_t remove_idx = dist(rng);
    int removed_city = new_path[remove_idx];

    //no longer visited
    visited[removed_city] = false;

    // Compute cost difference
    double cost_before = (remove_idx > 1) ? costs[new_path[remove_idx - 1]][removed_city] : 0.0;
    double cost_after = (remove_idx < new_path.size() - 1) ? costs[removed_city][new_path[remove_idx + 1]] : 0.0;
    double new_segment_cost = (remove_idx < new_path.size() - 1) ? costs[new_path[remove_idx - 1]][new_path[remove_idx + 1]] : 0.0;
    double freed_budget = cost_before + cost_after - new_segment_cost;

    new_path.erase(new_path.begin() + remove_idx);
    current_cost -= freed_budget;

    int current_city = new_path.back();
    int num_cities = static_cast<int>(costs.size());

    while (current_cost < budget) {
        int best_city = -1;
        double min_cost = std::numeric_limits<double>::infinity();

        // Pick a city not yet visited
        for (int i = 0; i < num_cities; ++i) {
            if (!visited[i] && i != current_city) {
                double extra_cost = costs[current_city][i];
                if (current_cost + extra_cost <= budget && extra_cost < min_cost) {
                    min_cost = extra_cost;
                    best_city = i;
                }
            }
        }

        if (best_city == -1) break;

        new_path.push_back(best_city);
        visited[best_city] = true;   //visited
        current_cost += min_cost;
        current_city = best_city;
    }

    return new_path;
}


std::vector<int> repair_path(const std::vector<int>& path,
                             const std::vector<std::vector<double>>& costs,
                             const std::vector<double>& prizes,
                             double budget,
                             std::vector<bool>& visited) {
    std::vector<int> new_path = path;
    double current_cost = path_cost(costs, path);
    int num_cities = static_cast<int>(costs.size());
    int start_city = path.front();

    //remove cities with least prize if exceeding budget
    while (current_cost > budget && new_path.size() > 1) {
        int min_prize_city = -1;
        double min_prize = std::numeric_limits<double>::max();
        for (size_t i = 0; i < new_path.size(); ++i) {
            int city = new_path[i];
            if (city != start_city && prizes[city] < min_prize) {
                min_prize = prizes[city];
                min_prize_city = i;
            }
        }

        if (min_prize_city == -1) break;

        int removed_city = new_path[min_prize_city];
        new_path.erase(new_path.begin() + min_prize_city);
        visited[removed_city] = false;

        current_cost = path_cost(costs, new_path);
    }

    //add high-prize cities if within budget
    while (true) {
        std::vector<std::pair<int, double>> candidates;
        int current_city = new_path.back();

        for (int i = 0; i < num_cities; ++i) {
            if (!visited[i] && i != current_city) {
                double extra_cost = costs[current_city][i];
                if (current_cost + extra_cost <= budget) {
                    candidates.emplace_back(i, prizes[i]);
                }
            }
        }

        if (candidates.empty()) break;

        auto best_candidate = std::max_element(candidates.begin(), candidates.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

        int new_city = best_candidate->first;
        new_path.push_back(new_city);
        visited[new_city] = true;
        current_cost += costs[current_city][new_city];
    }

    return new_path;
}



std::vector<int> neighbor(const std::vector<int>& path,
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double budget,
                          std::vector<bool>& visited) {
    std::uniform_real_distribution<double> method_selector(0.0, 1.0);
    double method = method_selector(rng);

    std::vector<int> new_path;
    if (method < 0.2) {
        new_path = swap_cities(path); 
    }
    else if (method < 0.4) {
        new_path = change_city(path, costs, prizes, budget, visited);
    }
    else if (method < 0.6) {
        new_path = shift_segment(path);
    }
    else if (method < 0.8) {
        new_path = shuffle_segment(path);
    }
    else {
        new_path = two_opt(path);
    }

    return repair_path(new_path, costs, prizes, budget, visited);
}


std::vector<int> simulated_annealing(const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& prizes,
                                     double budget, int start_city,
                                     std::vector<int> init_path,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate) {
    if (init_path.empty() || init_path[0] != start_city) {
        init_path = {start_city};
    }

    auto current_path = init_path;
    auto best_path = current_path;
    double best_prize = path_prize(prizes, best_path);

    std::vector<bool> visited(costs.size(), false);
    for (int city : current_path) {
        visited[city] = true;
    }

    double temp = initial_temp;

    for (int iter = 0; iter < iterations; iter++) {
        auto local_visited = visited;
        auto candidate_path = neighbor(current_path, costs, prizes, budget, local_visited);
        if (candidate_path[0] != start_city) {
            continue;
        }

        double candidate_cost = path_cost(costs, candidate_path);
        if (candidate_cost > budget) {
            continue;
        }

        double candidate_prize = path_prize(prizes, candidate_path);
        double delta = candidate_prize - path_prize(prizes, current_path);

        if (delta >= 0 || exp(delta / temp) > std::uniform_real_distribution<>(0.0, 1.0)(rng)) {
            current_path = candidate_path;
            visited = local_visited;
            if (candidate_prize > best_prize) {
                best_path = candidate_path;
                best_prize = candidate_prize;
            }
        }

        temp *= cooling_rate;
        if (temp < 1e-6) break;
    }

    return best_path;
}

std::vector<int> simulated_annealing_parallel(const std::vector<std::vector<double>>& costs,
                                              const std::vector<double>& prizes,
                                              double budget, int start_city,
                                              std::vector<std::vector<int>> init_paths,
                                              int iterations, double initial_temp,
                                              double cooling_rate, int num_threads) {
    std::vector<std::vector<int>> paths(num_threads);
    std::vector<double> prizes_found(num_threads, 0.0);

    std::random_device rd;
    unsigned int base_seed = rd();

    #pragma omp parallel for num_threads(num_threads)
    for (int thread = 0; thread < num_threads; thread++) {
        std::mt19937 thread_rng(base_seed + thread);
        double thread_temp = initial_temp * (0.9 + 0.2 * (thread / double(num_threads)));

        paths[thread] = simulated_annealing(
            costs, prizes, budget, start_city, init_paths[thread],
            iterations, thread_temp, cooling_rate
        );
        prizes_found[thread] = path_prize(prizes, paths[thread]);
    }

    auto best_it = std::max_element(prizes_found.begin(), prizes_found.end());
    int best_idx = std::distance(prizes_found.begin(), best_it);

    return paths[best_idx];
}

std::vector<int> simulated_annealing_optimization(const std::vector<std::vector<double>>& costs,
                                                 const std::vector<double>& prizes,
                                                 double budget, int start_city,
                                                 std::vector<int> init_path) {
    int iterations = 100000;
    double initial_temp = 100.0;
    double cooling_rate = 0.99995;
    int num_threads = 20;

    std::vector<std::vector<int>> init_paths(num_threads);

    if (init_path.empty()) {
        for (int i = 0; i < num_threads; ++i) {
            init_paths[i] = initial_solution(costs, budget, start_city);
        }
    } else {
        for (int i = 0; i < static_cast<int>(num_threads * 0.8); ++i) {
            init_paths[i] = init_path;
        }
        for (int i = static_cast<int>(num_threads * 0.8); i < num_threads; ++i) {
            init_paths[i] = initial_solution(costs, budget, start_city);
        }
    }

    return simulated_annealing_parallel(costs, prizes, budget, start_city, init_paths,
                                        iterations, initial_temp, cooling_rate, num_threads);
}