// #include "helplers.h"
#include "christofides.h"


#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <cmath>
#include <random>
#include <chrono>
#include <stack>
#include <unordered_set>
#include <stack>
#include <unordered_map>
#include <algorithm>

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


//Find odd-degree vertices
std::vector<int> findOddDegreeVertices(const std::list<std::shared_ptr<Edge>>& mst) {
    std::unordered_map<int, int> degree;
    for (const auto& edge : mst) {
        degree[edge->getHead()]++;
        degree[edge->getTail()]++;
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
            // matching.push_back(std::make_shared<Edge>(u, v, costs[u][v], 0));
            matching.push_back(std::make_shared<Edge>(u, v, costs[u][v]));

            used[i] = used[j] = true;
        }
    }

    return matching;
}


std::vector<std::shared_ptr<Edge>> combineGraph(const std::list<std::shared_ptr<Edge>>& mst, 
                                                const std::vector<std::shared_ptr<Edge>>& matching) {
    std::vector<std::shared_ptr<Edge>> eulerianGraph(mst.begin(), mst.end());
    eulerianGraph.insert(eulerianGraph.end(), matching.begin(), matching.end());
    return eulerianGraph;
}



std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start) {
    std::unordered_map<int, std::vector<int>> adjList;
    for (const auto& edge : eulerianGraph) {
        adjList[edge->getHead()].push_back(edge->getTail());
        adjList[edge->getTail()].push_back(edge->getHead());
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

// std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start) {
//     // Build adjacency list
//     std::unordered_map<int, std::vector<int>> adjList;
//     for (const auto& edge : eulerianGraph) {
//         adjList[edge->getHead()].push_back(edge->getTail());
//         adjList[edge->getTail()].push_back(edge->getHead());
//     }

//     // Check if start is valid
//     if (adjList.find(start) == adjList.end() || adjList[start].empty()) {
//         std::cerr << "Error: start vertex " << start << " is isolated or missing" << std::endl;
//         return {start};
//     }

//     std::vector<int> circuit;
//     std::stack<int> stack;
//     stack.push(start);

//     while (!stack.empty()) {
//         int curr = stack.top();
//         if (!adjList[curr].empty()) {
//             int next = adjList[curr].back();
//             adjList[curr].pop_back();
//             // Remove reverse edge
//             auto& revAdj = adjList[next];
//             revAdj.erase(std::find(revAdj.begin(), revAdj.end(), curr));
//             stack.push(next);
//         } else {
//             circuit.push_back(curr);
//             stack.pop();
//         }
//     }

//     // Reverse to get correct order
//     std::reverse(circuit.begin(), circuit.end());

//     // Verify circuit size
//     if (circuit.size() < 2) {
//         std::cerr << "Error: Eulerian circuit too small, size = " << circuit.size() << std::endl;
//     }

//     return circuit;
// }


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

void removeDummyEdges(std::list<std::shared_ptr<Edge>>& tree_edges, int num_real_nodes) {
    for (auto it = tree_edges.begin(); it != tree_edges.end(); ) {
        int h = (*it)->getHead();
        int t = (*it)->getTail();
        if (h >= num_real_nodes || t >= num_real_nodes) {
            it = tree_edges.erase(it);  // safe removal
        } else {
            ++it;
        }
    }
}


std::vector<int> christofides(const std::vector<std::vector<double>>& costs, 
                              std::list<std::shared_ptr<Edge>>& mst,
                              int start) {
    // Prune dummy edges
    removeDummyEdges(mst, costs.size());
    // std::cout << "num edges after prune: " << mst.size() << std::endl;
    // for(auto m : mst) {
    //     std:: cout << m->getHead() << "-" << m->getTail() << " ";
    // }
    // std::cout << std::endl;

    // // Verify start is in MST
    // std::unordered_set<int> mst_vertices;
    // for (const auto& edge : mst) {
    //     mst_vertices.insert(edge->getHead());
    //     mst_vertices.insert(edge->getTail());
    // }
    // if (mst_vertices.find(start) == mst_vertices.end()) {
    //     // Choose a valid start vertex (e.g., first vertex in MST)
    //     start = *mst_vertices.begin();
    //     std::cout << "Warning: Invalid start vertex, using " << start << std::endl;
    // }

    std::vector<int> oddVertices = findOddDegreeVertices(mst);
    // std::cout << "oddVertices: " << oddVertices.size() << std::endl;
    std::vector<std::shared_ptr<Edge>> matching = minimumWeightMatching(oddVertices, costs);
    // std::cout << "num matching: " << matching.size() << std::endl;
    std::vector<std::shared_ptr<Edge>> eulerianGraph = combineGraph(mst, matching);
    // std::cout << "num edges eulerianGraph: " << eulerianGraph.size() << std::endl;

    // // Verify start is in eulerianGraph
    // std::unordered_set<int> euler_vertices;
    // for (const auto& edge : eulerianGraph) {
    //     euler_vertices.insert(edge->getHead());
    //     euler_vertices.insert(edge->getTail());
    // }
    // if (euler_vertices.find(start) == euler_vertices.end()) {
    //     std::cerr << "Error: start vertex " << start << " not in eulerianGraph" << std::endl;
    //     return {};
    // }

    std::vector<int> circuit = findEulerianCircuit(eulerianGraph, start);
    // std::cout << "num vertices in circuit: " << circuit.size() << std::endl;
    // std::cout << "num edges in circuit: " << (circuit.size() > 1 ? circuit.size() - 1 : 0) << std::endl;
    std::vector<int> tspPath = shortcutHamiltonianPath(circuit);
    // std::cout << "tspPath vertices: " << tspPath.size() << std::endl;

    int last_pos = tspPath.back();
    if(costs[0][1] > costs[0][last_pos]) {
        std::reverse(tspPath.begin() + 1, tspPath.end());
    }
    return tspPath;
}
