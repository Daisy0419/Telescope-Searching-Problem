#pragma once

#include <memory>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <list>
#include <set>

class Edge {
public:
    int u, v;
    double weight;
    double potential;
    int active_ends;

    Edge(int u_, int v_, double weight_, int active_ends_);
    bool operator<=(const Edge& e) const;
};


class Subset {
public:
    int id;
    double dualValue;
    double potential;
    bool is_active;
    std::shared_ptr<Subset> root;
    
    std::unordered_set<int> vertices;
    std::pair<std::shared_ptr<Subset>, std::shared_ptr<Subset>> predecessors;
    std::shared_ptr<Edge> merge_edge;
    std::unordered_set<std::shared_ptr<Edge>> tied_edges;


    Subset(int vertex, double lambda);

    Subset(std::shared_ptr<Subset> s1, std::shared_ptr<Subset> s2, std::shared_ptr<Edge> e);

    bool operator<(const Subset& s) const;
};


struct CompareEdges {
    bool operator()(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs) const {
        double lhsRatio = (lhs->active_ends == 0) 
            ? std::numeric_limits<double>::infinity() 
            : (lhs->potential / lhs->active_ends);

        double rhsRatio = (rhs->active_ends == 0) 
            ? std::numeric_limits<double>::infinity() 
            : (rhs->potential / rhs->active_ends);

        const double EPS = 1e-12;
        if (std::fabs(lhsRatio - rhsRatio) > EPS) {
            return lhsRatio < rhsRatio;
        }
        
        return lhs.get() < rhs.get();
    }
};


bool isEdgeInList(const std::list<std::shared_ptr<Edge>>& edges, const std::shared_ptr<Edge>& edge);
void updateAndReorderEdges(std::set<std::shared_ptr<Edge>, CompareEdges>& active_edges, double delta_potential);
std::list<std::shared_ptr<Subset>> primalDualSubroutine(const std::vector<std::vector<double>>& costs, double lambda,
                                                     std::list<std::shared_ptr<Edge>>& forest, 
                                                     std::list<std::shared_ptr<Subset>>& tiedSubsets);
void prune(std::list<std::shared_ptr<Edge>>& forest, std::shared_ptr<Subset> subset);
void pruneAll(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets);
std::unordered_set<int> findLargestConnectedComponent(const std::list<std::shared_ptr<Edge>>& forest);
int revisitTiedEdges(std::list<std::shared_ptr<Edge>>& forest, std::vector<std::shared_ptr<Subset>>& tiedSubsets, 
                        std::list<std::shared_ptr<Subset>>& rootSubsets, size_t i);

// void revisitTiedEdges(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& tiedSubsets, 
//                         std::list<std::shared_ptr<Subset>>& rootSubsets);
std::shared_ptr<Subset> findSmallestSubsetContainingComponent(std::shared_ptr<Subset> s, 
                                        const std::unordered_set<int>& largestComponent);

int LargestKernel(std::list<std::shared_ptr<Edge>>& forest, std::shared_ptr<Subset> s);
int findLargestKernel(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets);
std::vector<std::shared_ptr<Edge>> buildMST(std::list<std::shared_ptr<Edge>>& forest, std::vector<int>& vertices);
bool hasIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2);
std::shared_ptr<Subset> pruneSunset(std::shared_ptr<Subset> s, const std::unordered_set<int>& largestComponent);
void pick(std::shared_ptr<Subset> s, int v, int k, std::vector<int>& pickedVertices);
std::vector<std::shared_ptr<Edge>> FindMST(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets, int k);
double findLambda1(const std::vector<std::vector<double>>& costs, int k, double epsilon = 1e-3);
// double findLambda1(const std::vector<std::vector<double>>& costs, std::list<std::shared_ptr<Edge>>& forest,
//                     std::list<std::shared_ptr<Subset>>& subsets, std::list<std::shared_ptr<Subset>>& tiedSubsets, 
//                     int k, double epsilon = 1e-3);
// void displayForest(const std::list<std::shared_ptr<Edge>>& forest);
std::vector<std::shared_ptr<Edge>> kMST(const std::vector<std::vector<double>>& costs, int k);
std::vector<int> findK(const std::vector<std::vector<double>>& costs, double budget);
void findBestTour(const std::vector<std::vector<double>>& costs, double budget);
