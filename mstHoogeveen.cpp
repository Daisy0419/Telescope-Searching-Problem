#include "mst.h"
#include "ReadData.h"
#include "helplers.h"

#include <bits/stdc++.h>
#include <future>

// Find vertices of "wrong degree" for one fixed endpoint s
std::vector<int> findWrongDegreeVertices(const std::list<std::pair<int, int>>& mst, int s, int numNodes) {
    std::unordered_map<int, int> degree;
    for (const auto& [u, v] : mst) {
        degree[u]++;
        degree[v]++;
    }
    std::vector<int> S;
    for (int v = 0; v < numNodes; ++v) {
        if (v == s) {
            if (degree[v] % 2 == 0) { // s has even degree
                S.push_back(v);
            }
        } else {
            if (degree[v] % 2 == 1) { // v != s has odd degree
                S.push_back(v);
            }
        }
    }
    return S;
}

// Greedy matching (leaves one vertex exposed if |S| is odd)
std::list<std::pair<int, int>> greedyMatching(const std::vector<int>& nodes,
                                                 const std::vector<std::vector<double>>& costMatrix) {
    std::list<std::pair<int, int>> matching;
    std::unordered_set<int> used;

    for (size_t i = 0; i < nodes.size(); ++i) {
        int u = nodes[i];
        if (used.count(u)) continue;

        double minCost = std::numeric_limits<double>::infinity();
        int bestV = -1;

        for (size_t j = i + 1; j < nodes.size(); ++j) {
            int v = nodes[j];
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
        }
    }
    return matching;
}

// Christofides heuristic for Hamiltonian path with one fixed endpoint s
std::vector<int> christofidesPathOneFixed(int numNodes,
                                          const std::list<std::pair<int, int>>& mstEdges,
                                          const std::vector<std::vector<double>>& costMatrix,
                                          int s) {
    auto S = findWrongDegreeVertices(mstEdges, s, numNodes);
    auto matchingEdges = greedyMatching(S, costMatrix);

    std::list<std::pair<int, int>> combined = mstEdges;
    combined.insert(combined.end(), matchingEdges.begin(), matchingEdges.end());

    auto adj = buildAdjacencyList(combined);

    //Ensure exactly two odd-degree vertices
    int oddCount = 0;
    for (const auto& [v, neighbors] : adj) {
        if (neighbors.size() % 2 == 1) {
            oddCount++;
        }
    }
    if (oddCount == 0) {
        // All degrees even: delete an edge incident to s
        if (!adj[s].empty()) {
            int w = adj[s].front();
            adj[s].remove(w);
            adj[w].remove(s);
        } else {
            throw std::runtime_error("Fixed endpoint s has no neighbors in the combined graph");
        }
    }

    auto eulerPath = findEulerianTour(adj, s);

    return shortcutTour(eulerPath, s);
}

void binarySearchBestPathHoogeveen(const std::vector<std::vector<double>>& costs, 
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

        std::vector<int> mst_path = christofidesPathOneFixed(N, mst, costs, indices[0]);

        fix_cross(mst_path, costs);

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path.size(); ++i) {
            pathCost += costs[mst_path[i-1]][mst_path[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path) {
                new_prize += prizes[node];
            }

            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path;
            }
            N_min = N;
        }
    }
}

std::vector<int> mstHoogeveen(const std::vector<std::vector<double>>& costs, 
                           const std::vector<double>& prizes, 
                           double budget, 
                           int init_pos)  {
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
    binarySearchBestPathHoogeveen(costs, prizes, indices, best_mst_path, N_min, N_max, budget);

    while (best_mst_path.empty() && N_min >= 2) {
        N_max = N_min;
        N_min /= 2;
        binarySearchBestPathHoogeveen(costs, prizes, indices, best_mst_path, N_min, N_max, budget);
    }

    return best_mst_path;
}


std::vector<int> christofidesPathOneFixed(const Graph& g, const std::vector<std::vector<double>>& costs,
                                             const std::vector<int>& indices, int N_nodes,
                                             std::vector<Edge>& mst_edges, int fixed_start_idx) {
    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    // Determine vertices of wrong degree
    std::vector<int> wrong_vertices;
    for (int i = 0; i < N_nodes; ++i) {
        bool is_fixed = (i == fixed_start_idx);
        if ((is_fixed && degree[i] % 2 == 0) || (!is_fixed && degree[i] % 2 != 0)) {
            wrong_vertices.push_back(i);
        }
    }

    // Construct matching graph with dummy node if needed
    Graph odd_g;
    std::vector<Node> odd_nodes;
    int m = wrong_vertices.size();
    bool needs_dummy = (m % 2 != 0);

    for (int i = 0; i < m; ++i)
        odd_nodes.push_back(odd_g.addNode());

    Node dummy_node;
    if (needs_dummy) {
        dummy_node = odd_g.addNode();
        odd_nodes.push_back(dummy_node);
    }

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -costs[indices[wrong_vertices[i]]][indices[wrong_vertices[j]]];
        }
        if (needs_dummy) {
            Edge e = odd_g.addEdge(odd_nodes[i], dummy_node);
            odd_weight[e] = -1e6;  // large negative weight
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Build multigraph (MST + matching edges)
    std::vector<std::vector<int>> multigraph(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e));
        int v = g.id(g.v(e));
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    }

    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
        if (matching.matching(e)) {
            int u_idx = odd_g.id(odd_g.u(e));
            int v_idx = odd_g.id(odd_g.v(e));
            // std::cout << "u_idx: " << u_idx << ", v_idx: " << v_idx << "\n";
            if (needs_dummy && (u_idx == m || v_idx == m)) continue; // ignore dummy edges
            int u = wrong_vertices[u_idx];
            int v = wrong_vertices[v_idx];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    // std::cout << "multigraph side: " << multigraph.size() << std::endl;
    // Find Eulerian path starting from fixed_start_idx
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;
    std::vector<int> stack = {fixed_start_idx};

    while (!stack.empty()) {
        int v = stack.back();
        if (!temp_graph[v].empty()) {
            int u = temp_graph[v].back();
            // if (u >= N_nodes || v >= N_nodes) {
            //     std::cerr << "ERROR: out-of-bounds node: u=" << u << ", v=" << v << std::endl;
            // }
            temp_graph[v].pop_back();
            auto it = std::find(temp_graph[u].begin(), temp_graph[u].end(), v);
            if (it != temp_graph[u].end()) temp_graph[u].erase(it);
            stack.push_back(u);
        } else {
            circuit.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(circuit.begin(), circuit.end());
    // std::cout << "circuit.size()" << circuit.size() << std::endl;

    // Shortcut
    std::vector<bool> visited(N_nodes, false);
    std::vector<int> path;
    for (int v : circuit) {
        if (!visited[v]) {
            path.push_back(v);
            visited[v] = true;
        }
    }
    // std::cout << "path.size()" << path.size() << std::endl;

    return path;
}

void binarySearchBestPathOneFixed(const std::vector<std::vector<double>>& costs, 
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
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        std::vector<int> mst_path = std::move(christofidesPathOneFixed(g, costs, indices, N, mst_edges, 0));
        std::vector<int> mst_path_convert;
        for (int node : mst_path) {
            mst_path_convert.push_back(indices[node]);
        }        
        fix_cross(mst_path_convert, costs);

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
            //update best path
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path_convert;
            }

            N_min = N;
        }
    }
}


//s- path
std::vector<int> mstLemonHoogeveen(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
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
    
    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


// std::vector<int> christofidesPathTwoFixed(const Graph& g,
//                                           const std::vector<std::vector<double>>& costs,
//                                           const std::vector<int>& indices, int N_nodes,
//                                           const std::vector<Edge>& mst_edges,
//                                           int s_idx, int t_idx) {

//     std::vector<int> degree(N_nodes, 0);
//     for (Edge e : mst_edges) {
//         degree[g.id(g.u(e))]++;
//         degree[g.id(g.v(e))]++;
//     }

//     //vertices of “wrong” degree 
//     std::vector<int> wrong_vertices;
//     for (int i = 0; i < N_nodes; ++i) {
//         bool is_endpoint = (i == s_idx || i == t_idx);
//         bool should_be_odd = is_endpoint;        // s,t must be odd in final graph
//         if ((should_be_odd && degree[i] % 2 == 0) ||   // endpoint but even
//             (!should_be_odd && degree[i] % 2 != 0))    // interior but odd
//         {
//             wrong_vertices.push_back(i);
//         }
//     }

//     //minimum‑weight perfect matching on wrong vertices
//     Graph odd_g;
//     std::vector<Node> odd_nodes;
//     for (size_t k = 0; k < wrong_vertices.size(); ++k)
//         odd_nodes.push_back(odd_g.addNode());

//     Graph::EdgeMap<double> odd_weight(odd_g);
//     for (size_t i = 0; i < wrong_vertices.size(); ++i) {
//         for (size_t j = i + 1; j < wrong_vertices.size(); ++j) {
//             Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
//             odd_weight[e] = -costs[indices[wrong_vertices[i]]]
//                                    [indices[wrong_vertices[j]]];  // negate
//         }
//     }
//     lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>>
//         matching(odd_g, odd_weight);
//     matching.run();

//     std::vector<std::vector<int>> multigraph(N_nodes);
//     auto add_multiedge = [&](int u, int v) {
//         multigraph[u].push_back(v);
//         multigraph[v].push_back(u);
//     };
//     for (Edge e : mst_edges)
//         add_multiedge(g.id(g.u(e)), g.id(g.v(e)));

//     for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e)
//         if (matching.matching(e))
//             add_multiedge(wrong_vertices[odd_g.id(odd_g.u(e))],
//                           wrong_vertices[odd_g.id(odd_g.v(e))]);

//     //eulerian path from s_idx to t_idx
//     std::vector<int> circuit, stack = {s_idx};
//     std::vector<std::vector<int>> tmp = multigraph;
//     while (!stack.empty()) {
//         int v = stack.back();
//         if (!tmp[v].empty()) {
//             int u = tmp[v].back();
//             tmp[v].pop_back();
//             auto it = std::find(tmp[u].begin(), tmp[u].end(), v);
//             if (it != tmp[u].end()) tmp[u].erase(it);
//             stack.push_back(u);
//         } else {
//             circuit.push_back(v);
//             stack.pop_back();
//         }
//     }
//     std::reverse(circuit.begin(), circuit.end());

//     std::vector<bool> seen(N_nodes, false);
//     std::vector<int> path;
//     for (int v : circuit)
//         if (!seen[v]) { path.push_back(v); seen[v] = true; }

//     if (path.front() != s_idx) std::reverse(path.begin(), path.end());
//     return path;
// }


// std::vector<int> christofidesPathTwoFixed(const Graph& g,
//                                           const std::vector<std::vector<double>>& costs,
//                                           const std::vector<int>& indices, int N_nodes,
//                                           const std::vector<Edge>& mst_edges,
//                                           int s_idx, int t_idx) {
//     // Step 1: Find path from s to t in T and remove longest edge
//     std::vector<std::vector<int>> adj(N_nodes);
//     for (Edge e : mst_edges) {
//         int u = g.id(g.u(e)), v = g.id(g.v(e));
//         adj[u].push_back(v);
//         adj[v].push_back(u);
//     }
//     std::vector<int> parent(N_nodes, -1);
//     std::queue<int> q;
//     q.push(s_idx);
//     while (!q.empty() && parent[t_idx] == -1) {
//         int u = q.front(); q.pop();
//         for (int v : adj[u]) {
//             if (parent[v] == -1 && v != s_idx) {
//                 parent[v] = u;
//                 q.push(v);
//             }
//         }
//     }
//     std::vector<std::pair<int, int>> path_edges;
//     int curr = t_idx;
//     while (curr != s_idx) {
//         int prev = parent[curr];
//         path_edges.emplace_back(prev, curr);
//         curr = prev;
//     }
//     double max_cost = -1;
//     std::pair<int, int> e_to_remove;
//     for (auto [u, v] : path_edges) {
//         double cost = costs[indices[u]][indices[v]];
//         if (cost > max_cost) {
//             max_cost = cost;
//             e_to_remove = {u, v};
//         }
//     }

//     // Step 2: Construct bar{T} = T ∪ {d} \ {e}
//     std::vector<std::pair<int, int>> barT_edges;
//     for (Edge e : mst_edges) {
//         int u = g.id(g.u(e)), v = g.id(g.v(e));
//         if (!(u == e_to_remove.first && v == e_to_remove.second) &&
//             !(u == e_to_remove.second && v == e_to_remove.first)) {
//             barT_edges.emplace_back(u, v);
//         }
//     }
//     barT_edges.emplace_back(s_idx, t_idx); // Add dummy edge d

//     // Step 3: Compute degrees in bar{T}
//     std::vector<int> degree(N_nodes, 0);
//     for (auto [u, v] : barT_edges) {
//         degree[u]++;
//         degree[v]++;
//     }

//     // Step 4: Odd-degree vertices in bar{T}
//     std::vector<int> odd_vertices;
//     for (int i = 0; i < N_nodes; ++i) {
//         if (degree[i] % 2 != 0) {
//             odd_vertices.push_back(i);
//         }
//     }

//     // Step 5: Perfect matching on odd vertices
//     Graph odd_g;
//     std::vector<Node> odd_nodes;
//     for (int v : odd_vertices) odd_nodes.push_back(odd_g.addNode());
//     Graph::EdgeMap<double> odd_weight(odd_g);
//     for (size_t i = 0; i < odd_vertices.size(); ++i) {
//         for (size_t j = i + 1; j < odd_vertices.size(); ++j) {
//             Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
//             odd_weight[e] = -costs[indices[odd_vertices[i]]][indices[odd_vertices[j]]];
//         }
//     }
//     lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
//     matching.run();

//     // Step 6: Form multigraph
//     std::vector<std::vector<int>> multigraph(N_nodes);
//     auto add_multiedge = [&](int u, int v) {
//         multigraph[u].push_back(v);
//         multigraph[v].push_back(u);
//     };
//     for (auto [u, v] : barT_edges) add_multiedge(u, v);
//     for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
//         if (matching.matching(e)) {
//             int u = odd_vertices[odd_g.id(odd_g.u(e))];
//             int v = odd_vertices[odd_g.id(odd_g.v(e))];
//             add_multiedge(u, v);
//         }
//     }

//     // Step 7: Find Eulerian cycle
//     std::vector<int> circuit;
//     std::vector<int> stack = {s_idx};
//     std::vector<std::vector<int>> tmp = multigraph;
//     while (!stack.empty()) {
//         int v = stack.back();
//         if (!tmp[v].empty()) {
//             int u = tmp[v].back();
//             tmp[v].pop_back();
//             auto it = std::find(tmp[u].begin(), tmp[u].end(), v);
//             if (it != tmp[u].end()) tmp[u].erase(it);
//             stack.push_back(u);
//         } else {
//             circuit.push_back(v);
//             stack.pop_back();
//         }
//     }
//     // Ensure cycle by returning to start (since all degrees even)
//     circuit.push_back(s_idx);

//     // Step 8: Shortcut and remove dummy edge
//     std::vector<int> hamiltonian_cycle;
//     std::vector<bool> seen(N_nodes, false);
//     for (int v : circuit) {
//         if (!seen[v]) {
//             hamiltonian_cycle.push_back(v);
//             seen[v] = true;
//         }
//     }
//     // Find dummy edge d and ensure correct direction
//     std::vector<int> adjusted_cycle;
//     int d_pos = -1;
//     for (size_t i = 0; i < hamiltonian_cycle.size() - 1; ++i) {
//         int u = hamiltonian_cycle[i], v = hamiltonian_cycle[i + 1];
//         if ((u == s_idx && v == t_idx) || (u == t_idx && v == s_idx)) {
//             d_pos = i;
//             break;
//         }
//     }
//     if (hamiltonian_cycle[d_pos] == s_idx && hamiltonian_cycle[d_pos + 1] == t_idx) {
//         // d from s to t, reverse to get t to s
//         std::reverse(hamiltonian_cycle.begin(), hamiltonian_cycle.end());
//         d_pos = 0; // After reversal, adjust position if needed
//         for (size_t i = 0; i < hamiltonian_cycle.size() - 1; ++i) {
//             if (hamiltonian_cycle[i] == t_idx && hamiltonian_cycle[i + 1] == s_idx) {
//                 d_pos = i;
//                 break;
//             }
//         }
//     }
//     // Remove d and form path from s to t
//     std::vector<int> path;
//     for (size_t i = d_pos + 1; i < hamiltonian_cycle.size(); ++i) {
//         path.push_back(hamiltonian_cycle[i]);
//     }
//     for (size_t i = 0; i < d_pos; ++i) {
//         path.push_back(hamiltonian_cycle[i]);
//     }
//     return path;
// }


std::vector<int> christofidesPathTwoFixed(const Graph& g,
                         const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& indices,
                         int N_nodes,
                         const std::vector<Edge>& mst_edges,
                         int s_idx,
                         int t_idx) {
    //heaviest edge on the s-t path in MST
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
            odd_w[e] = -costs[indices[odd[i]]][indices[odd[j]]];
        }
    }

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
        // std::cout << "path_rev.front(): " << path_rev.front() << " path_rev.back(): " << path_rev.back();
        std::reverse(path_rev.begin(), path_rev.end());
        if (path_rev.front() != s_idx || path_rev.back() != t_idx)
            throw std::logic_error("Christofides 2-fixed: endpoints wrong.");
    }
    return path_rev;
}

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
        // fix_cross(mst_path_convert, costs);
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


//s-t path
std::vector<int> Hoogeveen2(const std::vector<std::vector<double>>& costs, 
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

    
    // std::cout << "Total elements: " << prizes.size() << std::endl;
    // std::cout << "N after phrase 1: " << N << std::endl;
    
    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


//select nodes by prize ratio w.r.t closest node already selected
void selectNodesPrizeRatio_(const std::vector<std::vector<double>>& costs,
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

std::vector<int> mstLemonHoogeveen2(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[init_pos] = true;
    indices.push_back(init_pos);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != init_pos) {
            minDist[v] = costs[init_pos][v];
            double score = std::asin(std::pow(prizes[v], 4)) / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 1, N = 2;
    while(true) {
        int num_to_select = N - N_pre;
        selectNodesPrizeRatio_(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if(N == (int)prizes.size())
            break;

        N_pre = N;
        N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }
    }

    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
                
}

//s-t path
std::vector<int> Hoogeveen22(const std::vector<std::vector<double>>& costs, 
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
        selectNodesPrizeRatio_(costs, prizes, indices, selected, minDist, heap, num_to_select);

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


//select nodes by prize ratio w.r.t closest node already selected
void selectNodesPrizeRatio_2(const std::vector<std::vector<double>>& costs,
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



//s-t path
std::vector<int> Hoogeveen23(const std::vector<std::vector<double>>& costs, 
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
        selectNodesPrizeRatio_2(costs, prizes, indices, selected, minDist, heap, num_to_select);

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


double sumPrizes(const std::vector<double>& prizes, const std::vector<int>& path) {
    double total_prizes = 0.0;
    for (int node : path) {
        total_prizes += prizes[node];
    }
    return total_prizes;
}

std::vector<int> bestHoogeveen(const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, 
                                    double budget, int s_idx, int t_idx) {

    auto f1 = std::async(std::launch::async, Hoogeveen2,  costs, prizes, budget, s_idx, t_idx);
    auto f2 = std::async(std::launch::async, Hoogeveen22, costs, prizes, budget, s_idx, t_idx);
    auto f3 = std::async(std::launch::async, Hoogeveen23, costs, prizes, budget, s_idx, t_idx);

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
