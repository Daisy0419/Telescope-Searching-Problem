#include "helplers.h"
#include "Christofides.h"


#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>
#include <stack>
#include <unordered_set>


UnionFind::UnionFind(const std::vector<int>& vertices) {
    num_elements = vertices.size();
    num_components = vertices.size();

    for(auto v : vertices) {
        rank[v] = 1;
        parent[v] = v;

    }
}

int UnionFind::find(int v) {
    int root = v;
    while (parent[root] != root) {
        root = parent[root];
    }
    // Path compression
    while (v != root) {
        int next = parent[v];
        parent[v] = root;
        v = next;
    }
    return root;
}

void UnionFind::unite(int u, int v) {
    int root_u = find(u);
    int root_v = find(v);
    
    if (root_u != root_v) {
        if (rank[root_u] >= rank[root_v]) {
            parent[root_v] = root_u;
            rank[root_u] += rank[root_v];
            
        } 
        else if (rank[root_u] < rank[root_v]) {
            parent[root_u] = root_v;
            rank[root_v] += rank[root_u];

        } 
        else {
            return;
        }

        num_components -= 1;
    }
}

int UnionFind::getSize(int x) {
    return rank[find(x)];
}


// //Kruskal's algorithm
// std::vector<SimpleEdge> minimumSpanningTree(const std::vector<std::vector<double>>& costs, std::vector<int>& vertices) {
//     std::vector<SimpleEdge> edges, mst;
//     int N = vertices.size();
//     UnionFind uf(vertices);

//     for (size_t i = 0; i < vertices.size(); i++) {
//         for (size_t j = i + 1; j < vertices.size(); j++) {
//             int u = vertices[i], v= vertices[j];
//             if (costs[u][v] != std::numeric_limits<double>::infinity()) {
//                 edges.push_back({u, v, costs[u][v]});
//             }
//         }
//     }

//     sort(edges.begin(), edges.end());

//     for (const auto& edge : edges) {
//         if (uf.find(edge.u) != uf.find(edge.v)) {
//             mst.push_back(edge);
//             uf.unite(edge.u, edge.v);
//         }
//     }
    
//     return mst;
// }

// //Find odd-degree vertices
// std::vector<int> findOddDegreeVertices(const std::vector<SimpleEdge>& mst) {
//     std::unordered_map<int, int> degree;
//     for (const auto& edge : mst) {
//         degree[edge.u]++;
//         degree[edge.v]++;
//     }

//     std::vector<int> oddVertices;
//     for (auto v : degree) {
//         if (v.second % 2 != 0) {
//             oddVertices.push_back(v.first);
//         }
//     }

//     std::cout << "num odd degree vertices: " << oddVertices.size() << std::endl;

//     return oddVertices;
// }

// //Min matching
// std::vector<SimpleEdge> minimumWeightMatching(const std::vector<int>& oddVertices, const std::vector<std::vector<double>>& costs) {
//     std::vector<SimpleEdge> matching;
//     std::vector<bool> used(oddVertices.size(), false);

//     for (size_t i = 0; i < oddVertices.size(); i++) {
//         int u = oddVertices[i];
//         if (used[i]) continue;

//         double min_cost = std::numeric_limits<double>::infinity();
//         int v = -1;
//         int j = 0;
//         for (; j < oddVertices.size(); j++) {
//             int candidate = oddVertices[j];
//             if (!used[j] && candidate != u && costs[u][candidate] < min_cost) {
//                 min_cost = costs[u][candidate];
//                 v = candidate;
//             }
//         }

//         if (v != -1) {
//             matching.push_back({u, v, costs[u][v]});
//             used[i] = used[j] = true;
//         }
//     }

//     return matching;
// }

// std::vector<SimpleEdge> combineGraph(const std::vector<SimpleEdge>& mst, const std::vector<SimpleEdge>& matching) {
//     std::vector<SimpleEdge> eulerianGraph = mst;
//     eulerianGraph.insert(eulerianGraph.end(), matching.begin(), matching.end());
//     return eulerianGraph;
// }


// std::vector<int> findEulerianCircuit(const std::vector<SimpleEdge>& eulerianGraph, int start) {
//     std::unordered_map<int, std::vector<int>> adjList;
//     for (const auto& edge : eulerianGraph) {
//         adjList[edge.u].push_back(edge.v);
//         adjList[edge.v].push_back(edge.u);
//     }

//     std::vector<int> circuit;
//     std::stack<int> stack;
//     stack.push(start);

//     while (!stack.empty()) {
//         int curr = stack.top();
//         if (!adjList[curr].empty()) {
//             int next = adjList[curr].back();
//             adjList[curr].pop_back();
//             auto& revAdj = adjList[next];
//             revAdj.erase(std::find(revAdj.begin(), revAdj.end(), curr));
//             stack.push(next);
//         } else {
//             circuit.push_back(curr);
//             stack.pop();
//         }
//     }
//     // std::cout << "vertices in Eulerian Circuit: " << circuit.size() << std::endl;
//     return circuit;
// }

// std::vector<int> shortcutHamiltonianPath(const std::vector<int>& eulerianCircuit) {
//     std::unordered_set<int> visited;
//     std::vector<int> tspPath;

//     for (int v : eulerianCircuit) {
//         if (visited.insert(v).second) {
//             tspPath.push_back(v);
//         }
//     }
//     // std::cout << "vertices in Hamiltonian Path: " << tspPath.size() << std::endl;
//     return tspPath;
// }

// // Christofides Algorithm
// std::vector<int> christofides(const std::vector<std::vector<double>>& costs, std::vector<int>& vertices) {
//     std::vector<SimpleEdge> mst = minimumSpanningTree(costs, vertices);
//     std::vector<int> oddVertices = findOddDegreeVertices(mst);
//     std::vector<SimpleEdge> matching = minimumWeightMatching(oddVertices, costs);
//     std::vector<SimpleEdge> eulerianGraph = combineGraph(mst, matching);
//     std::vector<int> circuit = findEulerianCircuit(eulerianGraph, vertices[0]);
//     std::vector<int> tspPath = shortcutHamiltonianPath(circuit);

//     std::cout << "vertices in final path: " << tspPath.size() << std::endl;
//     std::cout << "Approximate TSP Path (Christofides): ";
//     for (int v : tspPath) {
//         std::cout << v << " ";
//     }
//     std::cout << std::endl;
//     return tspPath;
// }


//Find odd-degree vertices
std::vector<int> findOddDegreeVertices(const std::vector<std::shared_ptr<Edge>>& mst) {
    std::unordered_map<int, int> degree;
    for (const auto& edge : mst) {
        degree[edge->u]++;
        degree[edge->v]++;
    }

    std::vector<int> oddVertices;
    for (auto v : degree) {
        if (v.second % 2 != 0) {
            oddVertices.push_back(v.first);
        }
    }

    // std::cout << "num odd degree vertices: " << oddVertices.size() << std::endl;

    return oddVertices;
}


std::vector<std::shared_ptr<Edge>> minimumWeightMatching(const std::vector<int>& oddVertices, 
                                                         const std::vector<std::vector<double>>& costs) {
    std::vector<std::shared_ptr<Edge>> matching;
    std::vector<bool> used(oddVertices.size(), false);

    for (size_t i = 0; i < oddVertices.size(); i++) {
        int u = oddVertices[i];
        if (used[i]) continue;

        double min_cost = std::numeric_limits<double>::infinity();
        int v = -1;
        size_t j = 0;
        for (size_t k = 0; k < oddVertices.size(); k++) {
            int candidate = oddVertices[k];
            if (!used[k] && candidate != u && costs[u][candidate] < min_cost) {
                min_cost = costs[u][candidate];
                v = candidate;
                j = k; // Store the index for marking as used
            }
        }

        if (v != -1) {
            // Add the matching edge with default potential and active_ends
            matching.push_back(std::make_shared<Edge>(u, v, costs[u][v], 0));
            used[i] = used[j] = true;
        }
    }

    return matching;
}


std::vector<std::shared_ptr<Edge>> combineGraph(const std::vector<std::shared_ptr<Edge>>& mst, 
                                                const std::vector<std::shared_ptr<Edge>>& matching) {
    std::vector<std::shared_ptr<Edge>> eulerianGraph = mst;
    eulerianGraph.insert(eulerianGraph.end(), matching.begin(), matching.end());
    return eulerianGraph;
}



std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start) {
    std::unordered_map<int, std::vector<int>> adjList;
    for (const auto& edge : eulerianGraph) {
        adjList[edge->u].push_back(edge->v);
        adjList[edge->v].push_back(edge->u);
    }

    std::vector<int> circuit;
    std::stack<int> stack;
    stack.push(start);

    while (!stack.empty()) {
        int curr = stack.top();
        if (!adjList[curr].empty()) {
            int next = adjList[curr].back();
            adjList[curr].pop_back();
            auto& revAdj = adjList[next];
            revAdj.erase(std::find(revAdj.begin(), revAdj.end(), curr));
            stack.push(next);
        } else {
            circuit.push_back(curr);
            stack.pop();
        }
    }
    // std::cout << "vertices in Eulerian Circuit: " << circuit.size() << std::endl;
    return circuit;
}

std::vector<int> shortcutHamiltonianPath(const std::vector<int>& eulerianCircuit) {
    std::unordered_set<int> visited;
    std::vector<int> tspPath;

    for (int v : eulerianCircuit) {
        if (visited.insert(v).second) {
            tspPath.push_back(v);
        }
    }
    // std::cout << "vertices in Hamiltonian Path: " << tspPath.size() << std::endl;
    return tspPath;
}

// Christofides Algorithm
std::vector<int> christofides(const std::vector<std::vector<double>>& costs, 
                              const std::vector<std::shared_ptr<Edge>>& mst,
                              int start) {
    std::vector<int> oddVertices = findOddDegreeVertices(mst);
    std::vector<std::shared_ptr<Edge>> matching = minimumWeightMatching(oddVertices, costs);
    std::vector<std::shared_ptr<Edge>> eulerianGraph = combineGraph(mst, matching);
    std::vector<int> circuit = findEulerianCircuit(eulerianGraph, start);
    std::vector<int> tspPath = shortcutHamiltonianPath(circuit);

    // std::cout << "vertices in final path: " << tspPath.size() << std::endl;
    // std::cout << "Approximate TSP Path (Christofides): ";
    // for (int v : tspPath) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    return tspPath;
}



// std::vector<std::vector<double>> initialize_distances(int NUM_CITIES) {
//     std::vector<std::vector<double>> costs(NUM_CITIES, std::vector<double>(NUM_CITIES, std::numeric_limits<double>::infinity()));
//     for (int i = 0; i < NUM_CITIES; i++) {
//         for (int j = i + 1; j < NUM_CITIES; j++) {
//             costs[i][j] = random_double(10.0, 100.0); 
//             costs[j][i] =  costs[i][j];
//         }
//         costs[i][i] = std::numeric_limits<double>::infinity(); 
//     }
//     return costs;
// }


// int main(){
//     auto start = std::chrono::high_resolution_clock::now();
//     std::vector<std::vector<double>> costs = initialize_distances(200);
    
//     // double best_prize = 0;
//     // // double budget = 47358.8;

//     // double budget = 2000;
//     std::vector<int> vertices = unique_random_ints(0, 199, 50);
//     std::cout << "init vertices: " << vertices.size() << std::endl;

//     std::cout << " Initial Prize: " << vertices.size() << std::endl;
//     std::cout << "Initial Tour: ";
//     for (int city : vertices) {
//         std::cout << city << " ";
//     }
//     std::cout << std::endl;

//     double total_cost = 0.0;
//     for (size_t i = 0; i < vertices.size() - 1; i++) {
//         total_cost += costs[vertices[i]][vertices[i + 1]];
//         std::cout << total_cost << " ";
//     }
//     std::cout << std::endl;
//     std::cout << "Initial Cost: " << total_cost << std::endl;
//     std::vector<int> best_tour = christofides(costs, vertices);

//     std::cout << " Final Best Prize: " << best_tour.size() << std::endl;
//     std::cout << "Best Tour: ";
//     for (int city : best_tour) {
//         std::cout << city << " ";
//     }
//     std::cout << std::endl;

//     total_cost = 0.0;
//     for (size_t i = 0; i < best_tour.size() - 1; i++) {
//         total_cost += costs[best_tour[i]][best_tour[i + 1]];
//         std::cout << total_cost << " ";
//     }
//     std::cout << std::endl;
//     std::cout << "Total Cost: " << total_cost << std::endl;
    
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

//     return 0;
// }