#pragma once

#include "Graph.hpp"
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


std::vector<int> findOddDegreeVertices(const std::vector<std::shared_ptr<Edge>>& mst);
std::vector<std::shared_ptr<Edge>> minimumWeightMatching(const std::vector<int>& oddVertices, 
                                                         const std::vector<std::vector<double>>& costs);
std::vector<std::shared_ptr<Edge>> combineGraph(const std::list<std::shared_ptr<Edge>>& mst, 
                                                const std::vector<std::shared_ptr<Edge>>& matching);
std::vector<int> findEulerianCircuit(const std::vector<std::shared_ptr<Edge>>& eulerianGraph, int start);      
std::vector<int> shortcutHamiltonianPath(const std::vector<int>& eulerianCircuit);
std::vector<int> christofides(const std::vector<std::vector<double>>& costs, 
                              std::list<std::shared_ptr<Edge>>& mst,
                              int start);    

void removeDummyEdges(std::list<std::shared_ptr<Edge>>& tree_edges, int num_real_nodes);                                                                                   