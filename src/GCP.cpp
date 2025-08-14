#include "GCP.h"
#include "ReadData.h"
#include "helplers.h"

#include <vector>
#include <list>
#include <queue>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <functional>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <stack>
#include <fstream>
#include <sstream>
#include <future>
#include <cmath>


// sort nodes in a given range by prize
void sortNodesInRange(const std::vector<double>& prizes, std::vector<int>& indices, std::pair<int, int> range) {
    int n = prizes.size();
    if (range.second <= 0 || range.first < 0 || range.first >= range.second) {
        std::cerr << "Invalid inputs to sortNodesInRange\n";
        return;
    }
    if (range.second > n) range.second = n;

    //partial partition
    std::nth_element(indices.begin() + range.first, indices.begin() + range.second, indices.end(),
                     [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

    //sort within [range.first, range.second)
    std::sort(indices.begin() + range.first, indices.begin() + range.second,
              [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

}

void buildGraphUpdate(const std::vector<std::vector<double>>& costs, 
            const std::vector<int>& indices, int N_nodes,
            Graph& g, Graph::EdgeMap<double>& weight) {
    std::vector<Node> nodes(N_nodes);
    for (int i = 0; i < N_nodes; ++i)
        nodes[i] = g.addNode();
    for (int i = 0; i < N_nodes; ++i) {
        for (int j = i + 1; j < N_nodes; ++j) {
            Edge e = g.addEdge(nodes[i], nodes[j]);
            weight[e] = costs[indices[i]][indices[j]];
        }
    }
}


// Build minimun spanning tree
void buildMST(const Graph& g, const Graph::EdgeMap<double>& weight, std::vector<Edge>& mst_edges) {
    mst_edges.clear();
    lemon::kruskal(g, weight, std::back_inserter(mst_edges));
}


// Compute Edge weight of a minimun spanning tree
double computeMSTWeight(const Graph& g, const Graph::EdgeMap<double>& weight, 
                        const std::vector<Edge>& mst_edges) {
    double total_weight = 0.0;
    for (const Edge& e : mst_edges) {
        total_weight += weight[e];
    }
    return total_weight;
}


//  Find odd-degree vertices
std::vector<int> findOddDegreeVertices(const std::list<std::pair<int, int>>& mst) {
    std::unordered_map<int, int> degree;
    for (const auto& [u, v] : mst) {
        degree[u]++;
        degree[v]++;
    }
    std::vector<int> odd;
    for (const auto& [node, deg] : degree) {
        if (deg % 2 == 1)
            odd.push_back(node);
    }
    return odd;
}

// Adapted christofides algorithm for a path with two fixed ends (J.A. Hoogeveen)
std::vector<int> christofidesPathTwoFixed(const Graph& g,
                         const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& indices,
                         int N_nodes,
                         const std::vector<Edge>& mst_edges,
                         int s_idx,
                         int t_idx) {

    //build adjacency lists
    std::vector<std::vector<int>> adj(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e)), v = g.id(g.v(e));
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    std::vector<int> parent(N_nodes, -1);
    std::queue<int> Q;
    Q.push(s_idx);
    parent[s_idx] = s_idx;
    while (!Q.empty() && parent[t_idx] == -1) {
        int u = Q.front(); Q.pop();
        for (int v : adj[u]) {
            if (parent[v] == -1) {
                parent[v] = u;
                Q.push(v);
            }
        }
    }

    // find heaviest edge on the s-t path in MST
    double max_cost = -1;
    std::pair<int, int> e_to_remove;
    for (int v = t_idx; v != s_idx; v = parent[v]) {
        int u = parent[v];
        double w = costs[indices[u]][indices[v]];
        if (w > max_cost) {
            max_cost = w;
            e_to_remove = {u, v};
        }
    }

    //Build T_bar = (T − {e_max}) ∪ {s,t}
    std::vector<std::pair<int, int>> tree_edges;
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e)), v = g.id(g.v(e));
        if (!((u == e_to_remove.first && v == e_to_remove.second) ||
              (u == e_to_remove.second && v == e_to_remove.first)))
        {
            tree_edges.emplace_back(u, v);
        }
    }
    tree_edges.emplace_back(s_idx, t_idx);  // dummy edge

    //find odd-degree vertices
    std::vector<int> degree(N_nodes, 0);
    for (auto [u, v] : tree_edges) {
        degree[u]++;
        degree[v]++;
    }

    std::vector<int> odd;
    for (int i = 0; i < N_nodes; ++i)
        if (degree[i] % 2 != 0)
            odd.push_back(i);

    //compute minimum-weight perfect matching on odd vertices
    Graph odd_g;
    std::vector<Node> odd_nodes;
    for (int v : odd)
        odd_nodes.push_back(odd_g.addNode());

    Graph::EdgeMap<double> odd_w(odd_g);
    for (size_t i = 0; i < odd.size(); ++i) {
        for (size_t j = i + 1; j < odd.size(); ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_w[e] = -costs[indices[odd[i]]][indices[odd[j]]]; //negate weight
        }
    }

    //compute minimum-weight perfect matching (using negated weight)
    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_w);
    matching.run();

    //create multigraph with T_bar ∪ matching edges
    std::vector<std::vector<int>> multigraph(N_nodes);
    auto add = [&](int u, int v) {
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    };
    for (auto [u, v] : tree_edges) add(u, v);
    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e)
        if (matching.matching(e)) {
            int u = odd[odd_g.id(odd_g.u(e))];
            int v = odd[odd_g.id(odd_g.v(e))];
            add(u, v);
        }

    //compute eulerian cycle
    std::vector<int> cyc, stack = {s_idx};
    auto temp = multigraph;
    while (!stack.empty()) {
        int v = stack.back();
        if (!temp[v].empty()) {
            int u = temp[v].back();
            temp[v].pop_back();
            auto it = std::find(temp[u].begin(), temp[u].end(), v);
            if (it != temp[u].end()) temp[u].erase(it);
            stack.push_back(u);
        } else {
            cyc.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(cyc.begin(), cyc.end());
    cyc.pop_back(); // remove duplicate start

    //break Eulerian cycle at dummy edge (s,t)
    size_t cut = cyc.size();
    for (size_t i = 0; i < cyc.size(); ++i) {
        size_t j = (i + 1) % cyc.size();
        int u = cyc[i], v = cyc[j];
        if ((u == s_idx && v == t_idx) || (u == t_idx && v == s_idx)) {
            cut = i;
            break;
        }
    }
    if (cut == cyc.size()) throw std::logic_error("Dummy edge not found!");

    std::rotate(cyc.begin(), cyc.begin() + cut, cyc.end());

    //shortcut duplicate nodes (keep first occurrence)
    std::vector<bool> seen(N_nodes, false);
    std::vector<int> path_rev;
    for (int i = 0; i < (int)cyc.size(); ++i) {
        int v = cyc[i];
        if (!seen[v]) {
            path_rev.push_back(v);
            seen[v] = true;
        }
    }
    std::reverse(path_rev.begin() + 1, path_rev.end());


    if (path_rev.front() != s_idx || path_rev.back() != t_idx){
        std::reverse(path_rev.begin(), path_rev.end());
        if (path_rev.front() != s_idx || path_rev.back() != t_idx)
            throw std::logic_error("Christofides 2-fixed: endpoints wrong.");
    }
    return path_rev;
}


// Greedy extend the path by adding a node in the end before t
void greedyExtendPath(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes,
                      std::vector<int>& path, double& total_cost, double budget) {
    int end_city = path.back();
    std::vector<bool> visited(costs.size(), false);
    for (int city : path) visited[city] = true;

    while (true) {
        double max_prize = -1.0;
        int best_insert_pos = -1;
        int next_city = -1;

        // Always consider inserting before the final node
        int before_end = path[path.size() - 2];

        for (int i = 0; i < (int)costs.size(); ++i) {
            if (visited[i]) continue;

            double extra_cost = costs[before_end][i] + costs[i][end_city] - costs[before_end][end_city];
            if (total_cost + extra_cost <= budget && prizes[i] > max_prize) {
                max_prize = prizes[i];
                next_city = i;
            }
        }

        if (next_city == -1) break;

        // Insert next_city before end_city
        path.insert(path.end() - 1, next_city);
        visited[next_city] = true;

        // Recalculate last segment before end
        int before_end_new = path[path.size() - 2];
        int before_before = path[path.size() - 3];
        total_cost += costs[before_before][before_end_new] + costs[before_end_new][end_city]
                      - costs[before_before][end_city];
    }
}


void binarySearchBestPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    const std::vector<int>& indices,
                    std::vector<int>& best_mst_path,
                    int N_min, int N_max, double budget) {
    
    double best_prize = 0.0;
    N_max += 1;
    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        // std::cout << "N_max: " << N_max << " N_min: " << N_min << " N: " << N << std::endl;;
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        // std::cout << "MST weight = " << total_mst_weight << std::endl;
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        if (N == 2) {
            double path_cost = costs[indices[0]][indices[1]];
            if (path_cost <= budget) {
                best_mst_path = {indices[0], indices[1]};
            } else {
                best_mst_path = {indices[0]};
            }
            return;
        }
        std::vector<int> mst_path = std::move(christofidesPathTwoFixed(g, costs, indices, N, mst_edges, 0, 1));
        std::vector<int> mst_path_convert;
        for (int node : mst_path) {
            mst_path_convert.push_back(indices[node]);
        }        
        //2-opt to remove cross
        fix_cross_st_path(mst_path_convert, costs);


        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path_convert.size(); ++i) {
            pathCost += costs[mst_path_convert[i-1]][mst_path_convert[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            greedyExtendPath(costs, prizes, mst_path_convert, pathCost, budget);
            double new_prize = 0.0;
            for (int node : mst_path_convert) {
                new_prize += prizes[node];
            }
            // std::cout << "new_prize: " << new_prize << std::endl;
            //update best path
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path_convert;
                // std::cout << "best_prize = " << best_prize << std::endl;
            }

            N_min = N;
        }
    }
}


//s-t path - Best Prize Node selection with one fixed end Christofedis
std::vector<int> GCP_BestPrizeTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    std::vector<int> indices(prizes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::swap(indices[0], indices[s_idx]);
    std::swap(indices[1], indices[t_idx]);

    int N_pre = 2, N = 4;
    while(true) {
        sortNodesInRange(prizes, indices, {N_pre, N});

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N, N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }

    }

    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


// Select nodes using nonlinear prize boosting: score = asin(prize^4) / (minDist + epsilon)
// This emphasizes higher-prize nodes more aggressively and damps extreme values using asin().
void selectNodesByPrizeBoosting(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices, 
                            std::vector<bool>& selected,
                            std::vector<double>& minDist,
                            std::priority_queue<std::pair<double, int>>& heap,
                            int num) {

    int n = prizes.size();
    const double epsilon = 1;

    while (num > 0 && !heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);
        num --;
        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = std::asin(std::pow(prizes[v], 4))/ (minDist[v] + epsilon);
                heap.push({new_score, v});
            }
        }
    }

}

//s-t path
std::vector<int> GCP_PrizeBoostTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[s_idx] = true;
    indices.push_back(s_idx);
    selected[t_idx] = true;
    indices.push_back(t_idx);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != s_idx && v != t_idx ) {
            minDist[v] = std::min(costs[s_idx][v], costs[t_idx][v]);
            double score = std::asin(std::pow(prizes[v], 4)) / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 2, N = 4;
    while(true) {
        int num_to_select = N - N_pre;
        selectNodesByPrizeBoosting(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N, N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }

    }
    
    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


// Select nodes greedily based on prize-to-distance ratio: score = prize / (minDist + epsilon)
// Which favors nodes that are close to already selected nodes and have high prize values.
void selectNodesByPrizeDistanceRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices, 
                            std::vector<bool>& selected,
                            std::vector<double>& minDist,
                            std::priority_queue<std::pair<double, int>>& heap,
                            int num) {

    int n = prizes.size();
    const double epsilon = 1;

    while (num > 0 && !heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);
        num --;
        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = prizes[v]/ (minDist[v] + epsilon);
                heap.push({new_score, v});
            }
        } 
    }

}

// find best s-t path within budget
std::vector<int> GCP_PrizeRatioTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[s_idx] = true;
    indices.push_back(s_idx);
    selected[t_idx] = true;
    indices.push_back(t_idx);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != s_idx && v != t_idx ) {
            minDist[v] = std::min(costs[s_idx][v], costs[t_idx][v]);
            double score = prizes[v] / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 2, N = 4;
    while(true) {
        int num_to_select = N - N_pre;
        selectNodesByPrizeDistanceRatio(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N, N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }

    }
    
    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


// Computate total prize of a path
double sumPrizes(const std::vector<double>& prizes, const std::vector<int>& path) {
    double total_prizes = 0.0;
    for (int node : path) {
        total_prizes += prizes[node];
    }
    return total_prizes;
}


// GCP algorithm wrapper
std::vector<int> GCP(const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, 
                                    double budget, int s_idx, int t_idx) {

    auto f1 = std::async(std::launch::async, GCP_BestPrizeTwoFix,  costs, prizes, budget, s_idx, t_idx);
    auto f2 = std::async(std::launch::async, GCP_PrizeBoostTwoFix, costs, prizes, budget, s_idx, t_idx);
    auto f3 = std::async(std::launch::async, GCP_PrizeRatioTwoFix, costs, prizes, budget, s_idx, t_idx);

    std::vector<int> p1 = f1.get();
    std::vector<int> p2 = f2.get();
    std::vector<int> p3 = f3.get();

    double s1 = sumPrizes(prizes, p1);
    double s2 = sumPrizes(prizes, p2);
    double s3 = sumPrizes(prizes, p3);

    if (s2 > s1 && s2 > s3) {
        return p2;
    } else if (s3 > s1 && s3 > s2) {
        return p3;
    } else {
        return p1;
    }
}

