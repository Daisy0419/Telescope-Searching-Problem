#include "mst.h"
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

#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/matching.h>

//sort by prize
void sortNodes(const std::vector<double>& prizes, std::vector<double>& prizes_sorted, 
               std::vector<int>& indices, int init_pos) {
    int n = prizes.size();
    if (init_pos < 0 || init_pos >= n) {
        prizes_sorted.clear();
        indices.clear();
        return;
    }
    std::vector<int> indices_temp(n);
    std::iota(indices_temp.begin(), indices_temp.end(), 0); // 0-based

    std::sort(indices_temp.begin(), indices_temp.end(), 
              [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

    prizes_sorted.resize(n);
    indices.resize(n);
    prizes_sorted[0] = prizes[init_pos];
    indices[0] = init_pos;
    size_t j = 1;
    for (size_t i = 0; i < n; ++i) {
        if (indices_temp[i] == init_pos)
            continue;
        prizes_sorted[j] = prizes[indices_temp[i]];
        indices[j] = indices_temp[i];
        ++j;
    }
}


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

//select nodes by prize ratio w.r.t current position
void selectNodes(const std::vector<std::vector<double>>& costs, 
                            const std::vector<double>& prizes, 
                            std::vector<double>& prizes_selected, 
                            std::vector<std::vector<double>>& costs_selected, 
                            std::vector<int>& indices, int init_pos,
                            double cost_threshold) {
    int n = prizes.size();

    //if unrooted 
    if (init_pos == -1) {
        init_pos = std::max_element(prizes.begin(), prizes.end()) - prizes.begin();
    }

    std::vector<bool> selected(n, false);

    int current_pos = init_pos;
    indices.push_back(current_pos);
    prizes_selected.push_back(prizes[current_pos]);
    selected[current_pos] = true;

    while (true) {
        double best_score = -1.0;
        int best_tile = -1;

        for (int j = 0; j < n; ++j) {
            if (selected[j] || j == current_pos) continue;

            double cost = costs[current_pos][j];
            if (cost <= cost_threshold) {
                double score = prizes[j] / (cost + 1E-10);
                if (score > best_score) {
                    best_score = score;
                    best_tile = j;
                }
            }
        }

        if (best_tile == -1) break;

        cost_threshold -= costs[current_pos][best_tile];
        current_pos = best_tile;
        selected[current_pos] = true;
        indices.push_back(current_pos);
        prizes_selected.push_back(prizes[current_pos]);

        if (cost_threshold <= 0) break;
    }

    int num_nodes = indices.size();
    costs_selected.resize(num_nodes, std::vector<double>(num_nodes));
    for (int i = 0; i < num_nodes; ++i) {
        for (int j = 0; j < num_nodes; ++j) {
            costs_selected[i][j] = costs[indices[i]][indices[j]];
        }
    }
}

//select nodes by prize ratio w.r.t closest node already selected
void selectNodesPrizeRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<double>& prizes_selected, 
                            std::vector<std::vector<double>>& costs_selected, 
                            std::vector<int>& indices, 
                            int init_pos, double cost_threshold) {
    const int n = costs.size();
    const double INF = std::numeric_limits<double>::max();
    const double EPS = 1e-10;

    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, INF);

    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio

    selected[init_pos] = true;
    indices.push_back(init_pos);
    prizes_selected.push_back(prizes[init_pos]);

    for (int v = 0; v < n; ++v) {
        if (v != init_pos) {
            minDist[v] = costs[init_pos][v];
            double score = prizes[v] / (minDist[v] + EPS);
            heap.push({score, v});
        }
    }

    double total_cost = 0.0;
    while (total_cost < cost_threshold && !heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);
        prizes_selected.push_back(prizes[u]);
        total_cost += minDist[u];

        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = prizes[v] / (minDist[v] + EPS);
                heap.push({new_score, v});
            }
        }
    }
    int num_nodes = indices.size();
    costs_selected.resize(num_nodes, std::vector<double>(num_nodes));
    for (int i = 0; i < num_nodes; ++i) {
        for (int j = 0; j < num_nodes; ++j) {
            costs_selected[i][j] = costs[indices[i]][indices[j]];
        }
    }    
}


void chooseFrontNodes(const std::vector<int>& indices,
                    const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<double>>& costs_filtered, int topN) {
    if (topN < 0 || topN > (int)indices.size() || topN > (int)costs.size()) {
        costs_filtered.clear();
        return;
    }
    costs_filtered.clear();
    costs_filtered.resize(topN, std::vector<double>(topN));
    for (int i = 0; i < topN; ++i) {
        for (int j = 0; j < topN; ++j) {
            costs_filtered[i][j] = costs[indices[i]][indices[j]];
        }
    }
}
//overload
void chooseFrontNodes(const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<double>>& costs_filtered, int topN) {
    costs_filtered.clear();
    costs_filtered.resize(topN, std::vector<double>(topN));
    for (int i = 0; i < topN; ++i) {
        for (int j = 0; j < topN; ++j) {
            costs_filtered[i][j] = costs[i][j];
        }
    }
}

std::list<std::pair<int, int>> buildMSTPrim(const std::vector<std::vector<double>>& costs) {
    int n = costs.size();
    std::vector<bool> inMST(n, false);
    std::vector<double> minCost(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);

    minCost[0] = 0.0;
    using P = std::pair<double, int>;  // (cost, vertex)
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    pq.emplace(0.0, 0);

    std::list<std::pair<int, int>> mstEdges;

    while (!pq.empty()) {
        auto [cost, u] = pq.top(); pq.pop();
        if (inMST[u]) continue;

        inMST[u] = true;

        if (parent[u] != -1) {
            mstEdges.emplace_back(parent[u], u); 
        }

        for (int v = 0; v < n; ++v) {
            if (!inMST[v] && costs[u][v] < minCost[v]) {
                minCost[v] = costs[u][v];
                parent[v] = u;
                pq.emplace(minCost[v], v);
            }
        }
    }

    return mstEdges;
}

std::list<std::pair<int, int>> buildMSTKruskal(const std::vector<std::vector<double>>& costs) {
    int n = costs.size(); 
    std::vector<std::tuple<double, int, int>> edges; // (weight, u, v)
    std::vector<int> parent(n), rank(n, 0); // For Union-Find

    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }

    //all edges
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            edges.emplace_back(costs[i][j], i, j);
        }
    }

    std::sort(edges.begin(), edges.end());

    std::function<int(int)> find = [&](int x) -> int {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); 
        }
        return parent[x];
    };

    auto unite = [&](int x, int y) {
        int px = find(x), py = find(y);
        if (px == py) return false;
        if (rank[px] < rank[py]) {
            parent[px] = py;
        } else if (rank[px] > rank[py]) {
            parent[py] = px;
        } else {
            parent[py] = px;
            rank[px]++;
        }
        return true;
    };

    //build MST
    std::list<std::pair<int, int>> mst;
    for (const auto& [weight, u, v] : edges) {
        if (unite(u, v)) {
            mst.emplace_back(u, v);
        }
    }

    return mst;
}

std::list<std::pair<int, int>> buildMSTKruskalUpdate(const std::vector<std::vector<double>>& costs, 
                                                const std::vector<int>& indices, int N_nodes) {

    std::vector<std::tuple<double, int, int>> edges; // (weight, u, v)
    std::vector<int> parent(N_nodes), rank(N_nodes, 0); // For Union-Find

    for (int i = 0; i < N_nodes; ++i) {
        parent[i] = i;
    }

    //all edges
    for (int i = 0; i < N_nodes; ++i) {
        for (int j = i + 1; j < N_nodes; ++j) {
            edges.emplace_back(costs[indices[i]][indices[j]], i, j);
        }
    }

    std::sort(edges.begin(), edges.end());

    std::function<int(int)> find = [&](int x) -> int {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); 
        }
        return parent[x];
    };

    auto unite = [&](int x, int y) {
        int px = find(x), py = find(y);
        if (px == py) return false;
        if (rank[px] < rank[py]) {
            parent[px] = py;
        } else if (rank[px] > rank[py]) {
            parent[py] = px;
        } else {
            parent[py] = px;
            rank[px]++;
        }
        return true;
    };

    //build MST
    std::list<std::pair<int, int>> mst;
    for (const auto& [weight, u, v] : edges) {
        if (unite(u, v)) {
            mst.emplace_back(indices[u], indices[v]);
        }
    }

    return mst;
}

double mstCosts(const std::vector<std::vector<double>>& costs, 
                const std::list<std::pair<int, int>>& mst) {
    double totalCost = 0.0;
    for (const auto& [u, v] : mst) {
        totalCost += costs[u][v];
    }

    return totalCost;
}

std::unordered_map<int, std::list<int>> buildAdjacencyList(const std::list<std::pair<int, int>>& edges) {
    std::unordered_map<int, std::list<int>> adj;
    for (const auto& [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    return adj;
}

//find Eulerian tour
std::vector<int> findEulerianTour(std::unordered_map<int, std::list<int>>& adj, int start) {
    std::vector<int> tour;
    std::stack<int> stack;
    stack.push(start);

    while (!stack.empty()) {
        int v = stack.top();
        if (!adj[v].empty()) {
            int u = adj[v].front();
            stack.push(u);
            adj[v].remove(u); 
            adj[u].remove(v); 
        } else {
            tour.push_back(v);
            stack.pop();
        }
    }

    std::reverse(tour.begin(), tour.end()); 
    return tour;
}

//shortcut to form Hamiltonian circuit
std::vector<int> shortcutTour(const std::vector<int>& eulerTour, int startNode) {
    std::unordered_set<int> visited;
    std::vector<int> tspTour;

    for (int v : eulerTour) {
        if (!visited.count(v)) {
            tspTour.push_back(v);
            visited.insert(v);
        }
    }

    // tspTour.push_back(tspTour.front());
    return tspTour;
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

//Greedy min matching
std::list<std::pair<int, int>> greedyPerfectMatching(const std::vector<int>& oddNodes,
                                                     const std::vector<std::vector<double>>& costMatrix) {
    std::list<std::pair<int, int>> matching;
    std::unordered_set<int> used;

    for (size_t i = 0; i < oddNodes.size(); ++i) {
        int u = oddNodes[i];
        if (used.count(u)) continue;

        double minCost = std::numeric_limits<double>::infinity();
        int bestV = -1;

        for (size_t j = i + 1; j < oddNodes.size(); ++j) {
            int v = oddNodes[j];
            if (used.count(v)) continue;

            if (costMatrix[u][v] < minCost) {
                minCost = costMatrix[u][v];
                bestV = v;
            }
        }

        if (bestV != -1) {
            matching.emplace_back(u, bestV);
            used.insert(u);
            used.insert(bestV);
        } else {
            std::cerr << "Warning: no match found for odd-degree node " << u << std::endl;
        }
    }

    return matching;
}

//wrapper for naive christofides
std::vector<int> christofidesTSP(int numNodes,
                                  const std::list<std::pair<int, int>>& mstEdges,
                                  const std::vector<std::vector<double>>& costMatrix,
                                  int startNode = 0) {
    auto oddVertices = findOddDegreeVertices(mstEdges);
    auto matchingEdges = greedyPerfectMatching(oddVertices, costMatrix);

    std::list<std::pair<int, int>> combined = mstEdges;
    combined.insert(combined.end(), matchingEdges.begin(), matchingEdges.end());

    auto adj = buildAdjacencyList(combined);
    auto eulerTour = findEulerianTour(adj, startNode);
    return shortcutTour(eulerTour, startNode);
}

void binarySearchBestPathNaive(const std::vector<std::vector<double>>& costs, 
                               const std::vector<double>& prizes,
                               const std::vector<int>& indices,
                               std::vector<int>& best_mst_path,
                               int N_min, int N_max, double budget) {
    double best_prize = 0.0;

    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        std::list<std::pair<int, int>> mst = buildMSTKruskalUpdate(costs, indices, N);
        double mst_cost = mstCosts(costs, mst);

        if (mst_cost > budget) {
            N_max = N;
            continue;
        }

        std::vector<int> mst_path = christofidesTSP(N, mst, costs, indices[0]);
        std::vector<int> mst_path_convert = mst_path;
        // for (int node : mst_path) {
        //     mst_path_convert.push_back(indices[node]);
        // }

        fix_cross(mst_path_convert, costs);

        int start_pos = mst_path_convert[0];
        int second_pos = mst_path_convert[1];
        int last_pos   = mst_path_convert.back();

        double dist_second = costs[start_pos][second_pos];
        double dist_last   = costs[start_pos][last_pos];
        if (dist_second > dist_last) {
            std::reverse(mst_path_convert.begin() + 1, mst_path_convert.end());
        }

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path_convert.size(); ++i) {
            pathCost += costs[mst_path_convert[i-1]][mst_path_convert[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path_convert) {
                new_prize += prizes[node];
            }

            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path_convert;
            }

            N_min = N;
        }
    }
}

std::vector<int> mstNaiveUpdate(const std::vector<std::vector<double>>& costs, 
                           std::vector<double>& prizes, 
                           double budget, 
                           int init_pos) 
{
    if (init_pos < 0 || init_pos >= (int)prizes.size() || budget <= 0) {
        return {init_pos};
    }

    std::vector<int> indices(prizes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::swap(indices[0], indices[init_pos]);

    int N = 2;
    while (true) {
        sortNodesInRange(prizes, indices, {N/2, N});
        std::list<std::pair<int, int>> mst = buildMSTKruskalUpdate(costs, indices, N);
        double total_mst_weight = mstCosts(costs, mst);

        if (total_mst_weight > budget) {
            break;
        }

        N *= 2;
        if (N >= (int)prizes.size()) {
            N = prizes.size();
            break;
        }
    }

    int N_min = N / 4, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathNaive(costs, prizes, indices, best_mst_path, N_min, N_max, budget);

    while (best_mst_path.empty() && N_min >= 2) {
        N_max = N_min;
        N_min /= 2;
        binarySearchBestPathNaive(costs, prizes, indices, best_mst_path, N_min, N_max, budget);
    }

    return best_mst_path;
}

std::vector<int> mstNaive(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {
    if (init_pos < 0 || init_pos >= (int)prizes.size() || budget <= 0) {
        return {init_pos};
    }

    std::vector<double> prizes_sorted; 
    std::vector<int> indices;
    sortNodes(prizes, prizes_sorted, indices, init_pos); 

    int N_min = 1, N_max = prizes.size();
    std::vector<int> best_mst_path;
    double best_prize = 0.0;
    double epsilon = 0.1;

    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;

        std::vector<std::vector<double>> costs_N;
        chooseFrontNodes(indices, costs, costs_N, N);
        // std::cout << "costs_N shape: " << costs_N.size() << " " << costs_N[0].size() << std::endl;
        // std::list<std::pair<int, int>> mst = buildMSTPrim(costs_N);
        std::list<std::pair<int, int>> mst = buildMSTKruskal(costs_N);
        double mst_cost = mstCosts(costs_N, mst);
        if(mst_cost > budget) {
            N_max = N;
            continue;
        }

        std::vector<int> mst_path = christofidesTSP(N, mst, costs_N, 0);
        // std::vector<int> mst_path = mstToPath(mst, N, 0);
        fix_cross(mst_path, costs_N);
        int start_pos = mst_path[0];
        int second_pos = mst_path[1];
        int last_pos   = mst_path.back();

        double dist_second = costs_N[start_pos][second_pos];
        double dist_last   = costs_N[start_pos][last_pos];

        if (dist_second > dist_last)
            std::reverse(mst_path.begin() + 1, mst_path.end());

        if (mst_path.size() < 2) {
            N_max = N;
            continue;
        }

        // Compute path
        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path.size(); ++i) {
            pathCost += costs_N[mst_path[i-1]][mst_path[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path) {
                new_prize += prizes[indices[node]];
            }

            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path;
            }

            // //if budget is nearly fully utilized
            // if (budget - pathCost < epsilon) {
            //     break;
            // }

            N_min = N;
        }
    }

    // convert to original indices
    std::vector<int> best_path;
    for (int node : best_mst_path) {
        best_path.push_back(indices[node]);
    }

    return best_path;
}

typedef lemon::ListGraph Graph;
typedef Graph::Node Node;
typedef Graph::Edge Edge;

void buildGraph(const std::vector<std::vector<double>>& costs, Graph& g, Graph::EdgeMap<double>& weight) {
    int n = costs.size();
    std::vector<Node> nodes(n);
    for (int i = 0; i < n; ++i)
        nodes[i] = g.addNode();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Edge e = g.addEdge(nodes[i], nodes[j]);
            weight[e] = costs[i][j];
        }
    }
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

void buildMST(const Graph& g, const Graph::EdgeMap<double>& weight, std::vector<Edge>& mst_edges) {
    mst_edges.clear();
    lemon::kruskal(g, weight, std::back_inserter(mst_edges));
}

double computeMSTWeight(const Graph& g, const Graph::EdgeMap<double>& weight, 
                        const std::vector<Edge>& mst_edges) {
    double total_weight = 0.0;
    for (const Edge& e : mst_edges) {
        total_weight += weight[e];
    }
    return total_weight;
}


std::vector<int> christofidesPath(const Graph& g, const std::vector<std::vector<double>>& costs,
                                std::vector<Edge>& mst_edges) {
    int n = costs.size();

    //odd degree vertices
    std::vector<int> degree(n, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    std::vector<int> odd_vertices;
    for (int i = 0; i < n; ++i) {
        if (degree[i] % 2 != 0)
            odd_vertices.push_back(i);
    }

    //Minimum Weight Perfect Matching
    Graph odd_g;
    std::vector<Node> odd_nodes(odd_vertices.size());
    for (int i = 0; i < odd_vertices.size(); ++i)
        odd_nodes[i] = odd_g.addNode();

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < odd_vertices.size(); ++i) {
        for (int j = i + 1; j < odd_vertices.size(); ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -costs[odd_vertices[i]][odd_vertices[j]]; //negate for min matching
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Merge MST + Matching
    std::vector<std::vector<int>> multigraph(n);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e));
        int v = g.id(g.v(e));
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    }

    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
        if (matching.matching(e)) {
            int u = odd_vertices[odd_g.id(odd_g.u(e))];
            int v = odd_vertices[odd_g.id(odd_g.v(e))];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    // Find Eulerian tour
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;
    std::vector<int> stack;
    stack.push_back(0); // starting at node 0

    while (!stack.empty()) {
        int v = stack.back();
        if (!temp_graph[v].empty()) {
            int u = temp_graph[v].back();
            temp_graph[v].pop_back();
            auto it = std::find(temp_graph[u].begin(), temp_graph[u].end(), v);
            if (it != temp_graph[u].end())
                temp_graph[u].erase(it);
            stack.push_back(u);
        } else {
            circuit.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(circuit.begin(), circuit.end());

    //shortcut
    std::vector<bool> visited(n, false);
    std::vector<int> path;
    for (int v : circuit) {
        if (!visited[v]) {
            path.push_back(v);
            visited[v] = true;
        }
    }
    // path.push_back(path[0]); 
    return path;
}

std::vector<int> christofidesPathUpdate(const Graph& g, const std::vector<std::vector<double>>& costs,
                                const std::vector<int>& indices, int N_nodes,
                                std::vector<Edge>& mst_edges) {

    //odd degree vertices
    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    std::vector<int> odd_vertices;
    for (int i = 0; i < N_nodes; ++i) {
        if (degree[i] % 2 != 0)
            odd_vertices.push_back(i);
    }

    //Minimum Weight Perfect Matching
    Graph odd_g;
    std::vector<Node> odd_nodes(odd_vertices.size());
    for (int i = 0; i < odd_vertices.size(); ++i)
        odd_nodes[i] = odd_g.addNode();

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < odd_vertices.size(); ++i) {
        for (int j = i + 1; j < odd_vertices.size(); ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -costs[indices[odd_vertices[i]]][indices[odd_vertices[j]]]; //negate for min matching
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Merge MST + Matching
    std::vector<std::vector<int>> multigraph(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e));
        int v = g.id(g.v(e));
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    }

    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
        if (matching.matching(e)) {
            int u = odd_vertices[odd_g.id(odd_g.u(e))];
            int v = odd_vertices[odd_g.id(odd_g.v(e))];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    // Find Eulerian tour
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;
    std::vector<int> stack;
    stack.push_back(0); // starting at node 0

    while (!stack.empty()) {
        int v = stack.back();
        if (!temp_graph[v].empty()) {
            int u = temp_graph[v].back();
            temp_graph[v].pop_back();
            auto it = std::find(temp_graph[u].begin(), temp_graph[u].end(), v);
            if (it != temp_graph[u].end())
                temp_graph[u].erase(it);
            stack.push_back(u);
        } else {
            circuit.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(circuit.begin(), circuit.end());

    //shortcut
    std::vector<bool> visited(N_nodes, false);
    std::vector<int> path;
    for (int v : circuit) {
        if (!visited[v]) {
            path.push_back(v);
            visited[v] = true;
        }
    }
    // path.push_back(path[0]); 
    return path;
}

std::vector<int> mstLemon(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {
    if (init_pos < 0 || init_pos >= (int)prizes.size() || budget <= 0) {
        return {init_pos};
    }

    std::vector<double> prizes_sorted; 
    std::vector<int> indices;
    sortNodes(prizes, prizes_sorted, indices, init_pos); 
    // for(int i = 0; i < 50; ++i) {
    //     std::cout << prizes[indices[i]] << " ";
    // }
    // std::cout << std::endl;

    int N_min = 1, N_max = prizes_sorted.size();
    std::vector<int> best_mst_path;
    double best_prize = 0.0;
    double epsilon = 0.1;

    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        std::vector<std::vector<double>> costs_N;
        chooseFrontNodes(indices, costs, costs_N, N);
        // std::vector<int> mst_path = christofides(costs_N);
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraph(costs_N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        // std::cout << "MST weight = " << total_mst_weight << std::endl;
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        std::vector<int> mst_path = std::move(christofidesPath(g, costs_N, mst_edges));
        mst_path.pop_back();
        fix_cross(mst_path, costs_N);

        int start_pos = mst_path[0];
        int second_pos = mst_path[1];
        int last_pos   = mst_path.back();

        double dist_second = costs_N[start_pos][second_pos];
        double dist_last = costs_N[start_pos][last_pos];
        if (dist_second > dist_last)
            std::reverse(mst_path.begin() + 1, mst_path.end());
            

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path.size(); ++i) {
            pathCost += costs_N[mst_path[i-1]][mst_path[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path) {
                new_prize += prizes_sorted[node];
            }

            //update best path
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path;
            }

            // //stop if budget is nearly fully utilized
            // if (budget - pathCost < epsilon) {
            //     break;
            // }

            N_min = N;
        }
    }

    //convert to original indices
    std::vector<int> best_path;
    for (int node : best_mst_path) {
        best_path.push_back(indices[node]);
    }
    return best_path;
}

void binarySearchBestPath(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    const std::vector<int>& indices,
                    std::vector<int>& best_mst_path,
                    int N_min, int N_max, double budget) {
    
    double best_prize = 0.0;
    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
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
        std::vector<int> mst_path = std::move(christofidesPathUpdate(g, costs, indices, N, mst_edges));
        mst_path.pop_back();
        // std::cout << "mst path len = " << mst_path.size() << std::endl;
        std::vector<int> mst_path_convert;
        for (int node : mst_path) {
            mst_path_convert.push_back(indices[node]);
        }        
        fix_cross(mst_path_convert, costs);

        int start_pos = mst_path_convert[0];
        int second_pos = mst_path_convert[1];
        int last_pos   = mst_path_convert.back();

        double dist_second = costs[start_pos][second_pos];
        double dist_last = costs[start_pos][last_pos];
        if (dist_second > dist_last)
            std::reverse(mst_path_convert.begin() + 1, mst_path_convert.end());
            

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path_convert.size(); ++i) {
            pathCost += costs[mst_path_convert[i-1]][mst_path_convert[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
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

std::vector<int> mstLemonUpdate(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {
    if (init_pos < 0 || init_pos >= (int)prizes.size() || budget <= 0) {
        return {init_pos};
    }

    std::vector<int> indices(prizes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::swap(indices[0], indices[init_pos]);

    int N = 2;
    while(true) {
        sortNodesInRange(prizes, indices, {N/2, N});
        // for(int i = 0; i < N; ++i) {
        //     std::cout << prizes[indices[i]] << " ";
        // }
        // std::cout << std::endl;
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
            break;
        }
    }
    
    int N_min = N/4, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPath(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPath(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


std::vector<int> mstLemon2(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    double cost_threshold = budget * 3;
    selectNodes(costs, prizes, prizes_selected, costs_selected, indices, init_pos, cost_threshold);

    int N_min = 1, N_max = prizes_selected.size();
    std::vector<int> best_mst_path;
    double best_prize = 0.0;
    double epsilon = 0.1;

    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        std::vector<std::vector<double>> costs_N;
        chooseFrontNodes(costs_selected, costs_N, N);
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraph(costs_N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        // std::cout << "MST weight = " << total_mst_weight << std::endl;
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        std::vector<int> mst_path = std::move(christofidesPath(g, costs_N, mst_edges));
        // std::vector<int> mst_path = christofides(costs_N);
        mst_path.pop_back();
        fix_cross(mst_path, costs_N);//2-opt

        int start_pos = mst_path[0];
        int second_pos = mst_path[1];
        int last_pos   = mst_path.back();

        double dist_second = costs_N[start_pos][second_pos];
        double dist_last = costs_N[start_pos][last_pos];
        if (dist_second > dist_last)
            std::reverse(mst_path.begin() + 1, mst_path.end());
            

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path.size(); ++i) {
            pathCost += costs_N[mst_path[i-1]][mst_path[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path) {
                new_prize += prizes[indices[node]];
            }

            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path;
            }

            if (budget - pathCost < epsilon) {
                break;
            }

            N_min = N;
        }
    }

    //convert to original indices
    std::vector<int> best_path;
    for (int node : best_mst_path) {
        best_path.push_back(indices[node]);
    }

    return best_path;
                
}



std::vector<int> mstLemon3(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    double cost_threshold = budget * 2;
    selectNodesPrizeRatio(costs, prizes, prizes_selected, costs_selected, indices, init_pos, cost_threshold);
    

    int N_min = 1, N_max = prizes_selected.size();
    std::vector<int> best_mst_path;
    double best_prize = 0.0;
    double epsilon = 0.1;

    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        std::vector<std::vector<double>> costs_N;
        chooseFrontNodes(costs_selected, costs_N, N);
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraph(costs_N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        // std::cout << "MST weight = " << total_mst_weight << std::endl;
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        std::vector<int> mst_path = std::move(christofidesPath(g, costs_N, mst_edges));
        // std::vector<int> mst_path = christofides(costs_N);
        mst_path.pop_back();
        fix_cross(mst_path, costs_N);//2-opt

        int start_pos = mst_path[0];
        int second_pos = mst_path[1];
        int last_pos   = mst_path.back();

        double dist_second = costs_N[start_pos][second_pos];
        double dist_last = costs_N[start_pos][last_pos];
        if (dist_second > dist_last)
            std::reverse(mst_path.begin() + 1, mst_path.end());
            

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path.size(); ++i) {
            pathCost += costs_N[mst_path[i-1]][mst_path[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path) {
                new_prize += prizes_selected[node];
            }

            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path;
            }

            if (budget - pathCost < epsilon) {
                break;
            }
            N_min = N;
        }
    }

    // Convert best_mst_path to original indices
    std::vector<int> best_path;
    for (int node : best_mst_path) {
        best_path.push_back(indices[node]);
    }

    return best_path;
                
}
