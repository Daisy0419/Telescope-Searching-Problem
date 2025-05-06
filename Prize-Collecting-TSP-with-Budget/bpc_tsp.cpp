
#include "Graph.hpp"
#include "Subset.hpp"
#include "Subroutine.hpp"
#include "PD.hpp"
#include "ReadFile.hpp"
#include "christofides.h"
#include "bpc_tsp.h"

#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <queue>
#include <numeric>





// int main(int argc, char* argv[]){
//     // std::cout << "start running..." << std::endl;
//     auto start = std::chrono::high_resolution_clock::now();
    
//     std::string file = "../new_data.vrp";
    
//     Graph G;
//     double meanEdgeWeight;
//     int numNodes;
//     std::string name;
//     // read the graph from the file and get some stats from it
//     // graphFromFile(file,G,name,meanEdgeWeight,numNodes);
//     buildGraph(file, G, meanEdgeWeight, numNodes);
//     // std::cout<<"mean edge weight:"<<meanEdgeWeight<<"\n"<<"number of nodes:"<<numNodes<<std::endl;
//     std::cout<< "Graph: "<<G<<"endGraph"<<std::endl;
    
//     // double l = 0.0020055;
//     // double r = 0.0020057;
//     // std::list<std::shared_ptr<Subset>> subsetsL = growSubsets(G,l);
//     // std::list<std::shared_ptr<Subset>> subsetsR = growSubsets(G,r);
//     // double wplus = reverseDelete(subsetsL,false);
//     // double wminus = reverseDelete(subsetsR,false);
    
//     std::list<std::shared_ptr<Edge>> mst;
//     double mst_w = G.MST(mst);
//     double D = 0.5*mst_w;
//     std::cout << "D is " << D << "\n";
    
//     std::list<std::shared_ptr<Edge>> edges;
//     double upper = 0;
//     int recursions = 1;
//     double lambda;
//     bool found;
//     PD(G,D,edges,upper,recursions,lambda,found,true);
    
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

// }


// std::vector<int> BPC_TSP(std::string file, const std::vector<std::vector<double>> &costs, 
//             const std::vector<double>& prizes, double D, int init_pos){

//     auto start = std::chrono::high_resolution_clock::now();
    
//     Graph G;
//     double meanEdgeWeight;
//     int numNodes;
//     std::string name;
    
//     buildGraph(file, G, meanEdgeWeight, numNodes, init_pos); 
//     // std::cout<<"mean edge weight:"<<meanEdgeWeight<<"\n"<<"number of nodes:"<<numNodes<<std::endl;
//     std::cout<< "Graph: "<<G<<"endGraph"<<std::endl;
    
//     std::cout << "D is " << D << "\n";

//     int l = D, r = 2*D;
//     std::vector<int> best_path;
//     double best_prize = 0;
//     while( l < r) {
//         double budget = (l+r)/2;
//         std::list<std::shared_ptr<Edge>> edges;
//         double upper = 0;
//         int recursions = 1;
//         double lambda;
//         bool found;
//         PD(G,budget,edges,upper,recursions,lambda,found,true);
//         std::vector<int> path = christofides(costs, edges, init_pos);
//         double total_cost = costs[path.back()][path.front()];
//         double total_prize = path.back();
//         for (size_t i = 0; i < path.size() - 1; i++) {
//             total_cost += costs[path[i]][path[i + 1]];
//             total_prize += prizes[path[i]];
//         }

//         if(total_cost < D) {
//             if(total_prize > best_prize) {
//                 best_path = path;
//             }
//             l = int(budget);
//         }
//         else {
//             r = int(budget);
//         }
//     }

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

//     return best_path;
// }

// std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start) {
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
//         std::cerr << "Error: Eulerian circuit too small" << std::endl;
//     }

//     return circuit;
// }

void sortNodes_(const std::vector<double>& prizes, std::vector<double>& prizes_sorted, 
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
    // for (size_t i = 0; i < n; ++i) {
    //     if (indices_temp[i] == init_pos)
    //         continue;
    //     prizes_sorted[j] = prizes[indices_temp[i]];
    //     indices[j] = indices_temp[i];
    //     ++j;
    // }
    for (size_t i = 0; i < n && j < n; ++i) {
        if (indices_temp[i] == init_pos)
            continue;
        prizes_sorted[j] = prizes[indices_temp[i]];
        indices[j] = indices_temp[i];
        ++j;
    }
}

// void chooseTopNodes_(const std::vector<double>& prizes_sorted, const std::vector<int>& indices,
//                     const std::vector<std::vector<double>>& costs,
//                     std::vector<std::vector<double>>& costs_filtered, 
//                     std::vector<double>& prizes_filtered, int topN) {
//     if (topN < 0 || topN > (int)indices.size() || topN > (int)costs.size()) {
//         costs_filtered.clear();
//         return;
//     }
//     prizes_filtered.clear();
//     prizes_filtered.assign(prizes_sorted.begin(), prizes_sorted.begin() + topN);
//     costs_filtered.resize(topN, std::vector<double>(topN));
//     for (int i = 0; i < topN; ++i) {
//         for (int j = 0; j < topN; ++j) {
//             costs_filtered[i][j] = costs[indices[i]][indices[j]];
//         }
//     }
// }


void chooseTopNodes_(const std::vector<double>& prizes_sorted, const std::vector<int>& indices,
                     const std::vector<std::vector<double>>& costs,
                     std::vector<std::vector<double>>& costs_filtered, 
                     std::vector<double>& prizes_filtered, int topN) {
    if (topN < 0 || topN > (int)indices.size() || topN > (int)costs.size()) {
        costs_filtered.clear();
        prizes_filtered.clear();
        return;
    }

    // Validate indices
    for (int i = 0; i < topN; ++i) {
        if (indices[i] < 0 || indices[i] >= (int)costs.size()) {
            costs_filtered.clear();
            prizes_filtered.clear();
            throw std::out_of_range("Invalid index: " + std::to_string(indices[i]));
        }
    }

    // Clear and assign prizes_filtered
    prizes_filtered.clear();
    prizes_filtered.assign(prizes_sorted.begin(), prizes_sorted.begin() + topN);

    // Clear and resize costs_filtered to ensure topN × topN
    costs_filtered.clear();
    costs_filtered.resize(topN, std::vector<double>(topN));
    
    // Populate costs_filtered
    for (int i = 0; i < topN; ++i) {
        for (int j = 0; j < topN; ++j) {
            costs_filtered[i][j] = costs[indices[i]][indices[j]];
        }
    }

    // Verify costs_filtered is square
    for (const auto& row : costs_filtered) {
        if (row.size() != (size_t)topN) {
            throw std::runtime_error("costs_filtered row size mismatch: expected " + 
                                    std::to_string(topN) + ", got " + std::to_string(row.size()));
        }
    }
}



// std::list<std::pair<int, int>> buildMSTPrim_(const std::vector<std::vector<double>>& costs) {
//     int n = costs.size();
//     std::vector<bool> inMST(n, false);
//     std::vector<double> minCost(n, std::numeric_limits<double>::infinity());
//     std::vector<int> parent(n, -1);

//     minCost[0] = 0.0;
//     using P = std::pair<double, int>;  // (cost, vertex)
//     std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
//     pq.emplace(0.0, 0);

//     std::list<std::pair<int, int>> mstEdges;

//     while (!pq.empty()) {
//         auto [cost, u] = pq.top(); pq.pop();
//         if (inMST[u]) continue;

//         inMST[u] = true;

//         if (parent[u] != -1) {
//             mstEdges.emplace_back(parent[u], u);  // Add edge to MST
//         }

//         for (int v = 0; v < n; ++v) {
//             if (!inMST[v] && costs[u][v] < minCost[v]) {
//                 minCost[v] = costs[u][v];
//                 parent[v] = u;
//                 pq.emplace(minCost[v], v);
//             }
//         }
//     }

//     return mstEdges;
// }

std::list<std::pair<int, int>> buildMSTPrim_(const std::vector<std::vector<double>>& costs) {
    if (costs.empty()) {
        return {};
    }
    int n = costs.size();
    // Validate square matrix
    for (const auto& row : costs) {
        if (row.size() != n) {
            throw std::invalid_argument("Costs matrix is not square: expected " + 
                                       std::to_string(n) + ", got " + std::to_string(row.size()));
        }
    }

    std::vector<bool> inMST(n, false);
    std::vector<double> minCost(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);

    minCost[0] = 0.0;
    using P = std::pair<double, int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    pq.emplace(0.0, 0);

    std::list<std::pair<int, int>> mstEdges;

    while (!pq.empty()) {
        auto [cost, u] = pq.top(); pq.pop();
        if (u < 0 || u >= n) {
            throw std::out_of_range("Invalid vertex index: " + std::to_string(u));
        }
        if (inMST[u]) continue;

        inMST[u] = true;

        if (parent[u] != -1) {
            mstEdges.emplace_back(parent[u], u);
        }

        for (int v = 0; v < n; ++v) {
            if (!inMST[v]) {
                if (v < 0 || v >= n || u < 0 || u >= n) {
                    throw std::out_of_range("Invalid access: u=" + std::to_string(u) + ", v=" + std::to_string(v));
                }
                if (std::isnan(costs[u][v]) || std::isinf(costs[u][v])) {
                    continue; // Skip invalid costs
                }
                if (costs[u][v] < minCost[v]) {
                    minCost[v] = costs[u][v];
                    parent[v] = u;
                    pq.emplace(minCost[v], v);
                }
            }
        }
    }

    return mstEdges;
}


std::vector<int> BPC_TSP(const std::vector<std::vector<double>>& costs,
                         const std::vector<double>& prizes,
                         double D,
                         int init_pos) {

    auto start = std::chrono::high_resolution_clock::now();

    Graph G;
    double meanEdgeWeight;
    int numNodes;
    std::string name;

    // buildGraph(file, G, meanEdgeWeight, numNodes, init_pos);
    buildGraph(costs, prizes, G, meanEdgeWeight, numNodes, init_pos);
    std::cout << "Graph: " << G << "endGraph" << std::endl;
    std::cout << "D is " << D << "\n";

    // Binary search bounds
    int l = D, r = 2 * D;
    std::vector<int> best_path;
    double best_prize = 0;

    while (r - l > 1) {
        double budget = int((l + r) / 2.0);
        std::cout << "budget: " << budget << std::endl;

        std::list<std::shared_ptr<Edge>> edges;
        double upper = 0;
        int recursions = 1;
        double lambda;
        bool found;

        PD(G, budget, edges, upper, recursions, lambda, found, true);
        std::cout << "num edges: " << edges.size() << std::endl;

        std::vector<int> path = christofides(costs, edges, init_pos);

        if (path.empty()) {
            r = int(budget);
            continue;
        }

        // for(auto p : path) {
        //     std::cout << p << " ";
        // }
        // std::cout << std::endl;

        double total_cost = 0;
        for (size_t i = 0; i < path.size() - 1; i++) {
            total_cost += costs[path[i]][path[i + 1]];
        }

        total_cost += costs[path.back()][path.front()];

        double total_prize = 0;
        for (size_t i = 0; i < path.size(); i++) {
            total_prize += prizes[path[i]];
        }

        if (total_cost < D) {
            if (total_prize > best_prize) {
                best_prize = total_prize;
                best_path = path;
            }
            l = int(budget);
        } else {
            r = int(budget);
        }
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

    best_path.push_back(init_pos);
    return best_path;
}

double mstCosts_(const std::vector<std::vector<double>>& costs, 
                const std::list<std::pair<int, int>>& mst) {
    // std::list<std::pair<int, int>> mst = buildMSTPrim(costs);

    double totalCost = 0.0;
    for (const auto& [u, v] : mst) {
        totalCost += costs[u][v];
    }

    return totalCost;
}

// std::vector<int> BPC_TSP2(const std::vector<std::vector<double>>& costs,
//                          const std::vector<double>& prizes,
//                          double D,
//                          int init_pos) {

//     std::vector<double> prizes_sorted; 
//     std::vector<int> indices;
//     sortNodes_(prizes, prizes_sorted, indices, init_pos); 

//     int N_min = 1, N_max = prizes.size(), N;
//     std::vector<std::vector<double>> costs_N;
//     std::vector<double> prizes_N;
//     while (N_max - N_min > 1) {
//         N = (N_max + N_min) / 2;

//         chooseTopNodes_(prizes_sorted, indices, costs, costs_N, prizes_N, N);
//         std::list<std::pair<int, int>> mst = buildMSTPrim_(costs_N);

//         double cost = mstCosts_(costs_N, mst);

//         if(cost > 2*D) {
//             N_max = N;
//         }
//         else {
//             N_min = N;
//         }  
//     }

//     chooseTopNodes_(prizes_sorted, indices, costs, costs_N, prizes_N, N);

//     Graph G;
//     double meanEdgeWeight;
//     int numNodes;
//     std::string name;

//     buildGraph(costs_N, prizes_N, G, meanEdgeWeight, numNodes, init_pos);
//     std::cout << "Graph: " << G << "endGraph" << std::endl;
//     std::cout << "D is " << D << "\n";

//     // Binary search bounds
//     int l = D, r = 2 * D;
//     std::vector<int> best_path;
//     double best_prize = 0;

//     while (r - l > 1) {
//         double budget = int((l + r) / 2.0);
//         std::cout << "budget: " << budget << std::endl;

//         std::list<std::shared_ptr<Edge>> edges;
//         double upper = 0;
//         int recursions = 1;
//         double lambda;
//         bool found;

//         PD(G, budget, edges, upper, recursions, lambda, found, true);
//         std::cout << "num edges: " << edges.size() << std::endl;

//         std::vector<int> path = christofides(costs_N, edges, init_pos);

//         if (path.empty()) {
//             r = int(budget);
//             continue;
//         }

//         // for(auto p : path) {
//         //     std::cout << p << " ";
//         // }
//         // std::cout << std::endl;

//         double total_cost = 0;
//         for (size_t i = 0; i < path.size() - 1; i++) {
//             total_cost += costs[path[i]][path[i + 1]];
//         }

//         total_cost += costs[path.back()][path.front()];

//         double total_prize = 0;
//         for (size_t i = 0; i < path.size(); i++) {
//             total_prize += prizes[path[i]];
//         }

//         if (total_cost < D) {
//             if (total_prize > best_prize) {
//                 best_prize = total_prize;
//                 best_path = path;
//             }
//             l = int(budget);
//         } else {
//             r = int(budget);
//         }
//     }

//     std::vector<int> converted_best_path;
//     for (int node : best_path) {
//         converted_best_path.push_back(indices[node]);
//     }

//     return converted_best_path;
// }


std::vector<int> BPC_TSP2(const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double D,
                          int init_pos) {
    std::vector<double> prizes_sorted; 
    std::vector<int> indices;
    sortNodes_(prizes, prizes_sorted, indices, init_pos); 

    int N_min = 1, N_max = prizes.size(), N;
    std::vector<std::vector<double>> costs_N;
    std::vector<double> prizes_N;
    while (N_max - N_min > 2) {
        N = (N_max + N_min) / 2;
        chooseTopNodes_(prizes_sorted, indices, costs, costs_N, prizes_N, N);
        std::cout << "N: " << N << std::endl;
        std::list<std::pair<int, int>> mst = buildMSTPrim_(costs_N);
        std::cout << "mst size: " << mst.size() << std::endl;
        double cost = mstCosts_(costs_N, mst);
        std::cout << "cost: " << cost << std::endl;
        if (cost > 2 * D) {
            N_max = N;
        } else {
            N_min = N;
        }  
    }

    std::cout << "N: " << N_min << std::endl;
    chooseTopNodes_(prizes_sorted, indices, costs, costs_N, prizes_N, N_min);

    Graph G;
    double meanEdgeWeight;
    int numNodes;
    std::string name;

    buildGraph(costs_N, prizes_N, G, meanEdgeWeight, numNodes, 0);
    std::cout << "Graph: " << G << "endGraph" << std::endl;
    std::cout << "D is " << D << "\n";

    int l = D, r = 2 * D;
    std::vector<int> best_path;
    double best_prize = 0;

    while (r - l > 1) {
        double budget = int((l + r) / 2.0);
        std::cout << "budget: " << budget << std::endl;

        std::list<std::shared_ptr<Edge>> edges;
        double upper = 0;
        int recursions = 1;
        double lambda;
        bool found;

        PD(G, budget, edges, upper, recursions, lambda, found, true);
        std::cout << "num edges: " << edges.size() << std::endl;

        std::vector<int> path = christofides(costs_N, edges, 0);

        if (path.empty()) {
            r = int(budget);
            continue;
        }

        double total_cost = 0;
        for (size_t i = 0; i < path.size() - 1; i++) {
            total_cost += costs_N[path[i]][path[i + 1]];
        }
        total_cost += costs_N[path.back()][path.front()];

        double total_prize = 0;
        for (size_t i = 0; i < path.size(); i++) {
            total_prize += prizes_N[path[i]];
        }

        if (total_cost <= D) {
            if (total_prize > best_prize) {
                best_prize = total_prize;
                best_path = path;
            }
            l = int(budget);
        } else {
            r = int(budget);
        }
    }

    if (best_path.empty()) {
        return {init_pos};
    }
    std::vector<int> converted_best_path;
    for (int node : best_path) {
        converted_best_path.push_back(indices[node]);
    }
    return converted_best_path;
}