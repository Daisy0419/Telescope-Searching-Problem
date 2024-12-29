#pragma once
#include "kMST.h"

#include <unordered_map>
#include <vector>

class UnionFind {
private:
    std::unordered_map<int, int> parent;
    std::unordered_map<int, int> rank;

public:
    int num_elements;
    int num_components;

    UnionFind(const std::vector<int>& vertices);

    int find(int v);
    void unite(int u, int v);
    int getSize(int x);

};

// struct SimpleEdge {
//     int u, v;
//     double weight;
//     bool operator<(const SimpleEdge& e) const {
//         return weight < e.weight;
//     }
// };


// std::vector<SimpleEdge> minimumSpanningTree(const std::vector<std::vector<double>>& costs, std::vector<int>& vertices);
// std::vector<int> findOddDegreeVertices(const std::vector<SimpleEdge>& mst);
// std::vector<SimpleEdge> minimumWeightMatching(const std::vector<int>& oddVertices, const std::vector<std::vector<double>>& costs);
// std::vector<SimpleEdge> combineGraph(const std::vector<SimpleEdge>& mst, const std::vector<SimpleEdge>& matching);
// std::vector<int> findEulerianCircuit(const std::vector<SimpleEdge>& eulerianGraph, int start);
// std::vector<int> shortcutHamiltonianPath(const std::vector<int>& eulerianCircuit);
// std::vector<int> christofides(const std::vector<std::vector<double>>& costs, std::vector<int>& vertices);


std::vector<int> findOddDegreeVertices(const std::vector<std::shared_ptr<Edge>>& mst);
std::vector<std::shared_ptr<Edge>> minimumWeightMatching(const std::vector<int>& oddVertices, 
                                                         const std::vector<std::vector<double>>& costs);
std::vector<std::shared_ptr<Edge>> combineGraph(const std::vector<std::shared_ptr<Edge>>& mst, 
                                                const std::vector<std::shared_ptr<Edge>>& matching);
std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start);      
std::vector<int> shortcutHamiltonianPath(const std::vector<int>& eulerianCircuit);
std::vector<int> christofides(const std::vector<std::vector<double>>& costs, 
                              const std::vector<std::shared_ptr<Edge>>& mst,
                              int start);                                                                                        