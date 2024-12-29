#include "helplers.h"
#include "kMST.h"
#include "Christofides.h"

#include <iostream>
#include <vector>
#include <unordered_set>
#include <memory>
#include <list>
#include <algorithm>
#include <map>
#include <limits>
#include <set>
#include <random>
#include <deque>
#include <fstream>
#include <sstream>
#include <chrono>
#include <queue>
#include <iomanip>

Subset::Subset(int vertex, double lambda) {
    id = vertex;
    vertices.insert(vertex);
    is_active = true;
    potential = lambda;
    dualValue = 0;
    predecessors = {nullptr, nullptr};
    merge_edge = nullptr;
}

Subset::Subset(std::shared_ptr<Subset> s1, std::shared_ptr<Subset> s2, std::shared_ptr<Edge> e) {
    id = s1->id;
    vertices.insert(s1->vertices.begin(), s1->vertices.end());
    vertices.insert(s2->vertices.begin(), s2->vertices.end());
    potential = s1->potential + s2->potential;
    dualValue = 0;
    is_active = true;
    predecessors = {s1, s2};
    merge_edge = e;

}

bool Subset::operator<(const Subset& s) const {
    return potential < s.potential;
}


Edge::Edge(int u_, int v_, double weight_, int active_ends_) : u(u_), v(v_), weight(weight_), 
    potential(weight_), active_ends(active_ends_) {}

bool Edge::operator<=(const Edge& e) const {
    if(active_ends == 0) {
        return false;
    }
    else if (e.active_ends == 0) {
        return true;
    }
    return potential / active_ends < e.potential / e.active_ends;
}



bool isEdgeInList(const std::list<std::shared_ptr<Edge>>& edges, const std::shared_ptr<Edge>& edge) {
    return std::find(edges.begin(), edges.end(), edge) != edges.end();
}

void updateAndReorderEdges(std::set<std::shared_ptr<Edge>, CompareEdges>& active_edges, double delta_potential) {
    std::vector<std::shared_ptr<Edge>> temp_edges(active_edges.begin(), active_edges.end());
    active_edges.clear(); 

    for (auto& edge : temp_edges) {
        edge->potential -= delta_potential * edge->active_ends;
    }

    active_edges.insert(temp_edges.begin(), temp_edges.end());
}





std::list<std::shared_ptr<Subset>> primalDualSubroutine(const std::vector<std::vector<double>>& costs, double lambda,
                                                     std::list<std::shared_ptr<Edge>>& forest, 
                                                     std::list<std::shared_ptr<Subset>>& tiedSubsets) {
    int num_vertices = costs.size();
    std::map<int, std::shared_ptr<Subset>> vertexToSubset;
    std::vector<std::vector<std::shared_ptr<Edge>>> edges(num_vertices, std::vector<std::shared_ptr<Edge>>(num_vertices, nullptr));
    std::list<std::shared_ptr<Subset>> rootSubsets;

    // Initialize subsets, one for each vertex
    std::list<std::shared_ptr<Subset>> activeSubsets;
    for (int i = 0; i < num_vertices; ++i) {
        auto subset = std::make_shared<Subset>(i, lambda);
        activeSubsets.push_back(subset);
        rootSubsets.push_back(subset);
        vertexToSubset[i] = subset;
    }
    // std::cout << "initial num vertices: " << activeSubsets.size() << std::endl;

    // Initialize edges
    std::set<std::shared_ptr<Edge>, CompareEdges> active_edges;
    for (size_t i = 0; i < num_vertices; i++) {
        for (size_t j = i + 1; j < num_vertices; j++) {
            if (costs[i][j] != std::numeric_limits<double>::infinity()) {
                auto edge = std::make_shared<Edge>(i, j, costs[i][j], 2);
                active_edges.insert(edge);
                edges[i][j] = edge;
                edges[j][i] = edge;
            }
        }
    }
    // std::cout << "initial num edges: " << active_edges.size() << std::endl;

    while (activeSubsets.size() > 0) {
        std::shared_ptr<Edge> next_edge = nullptr;
    
        if (!active_edges.empty()) {
            next_edge = *active_edges.begin();
        }
        auto next_subset = activeSubsets.front();

        double delta_potential;
        if (next_edge != nullptr && next_edge->potential / next_edge->active_ends <= next_subset->potential) {

            delta_potential = next_edge->potential / next_edge->active_ends;
            for (auto& subset : activeSubsets) {
                subset->potential -= delta_potential;
            }

            updateAndReorderEdges(active_edges, delta_potential);

            forest.push_back(next_edge);
            active_edges.erase(next_edge);

            // handle subset
            auto s1 = vertexToSubset[next_edge->u];
            auto s2 = vertexToSubset[next_edge->v];

            auto newSubset = std::make_shared<Subset>(s1, s2, next_edge);
            activeSubsets.remove(s1);
            activeSubsets.remove(s2);

            for(int v : newSubset->vertices) {
                vertexToSubset[v] = newSubset;
            }

            auto it = std::find_if(activeSubsets.begin(), activeSubsets.end(), [&](const std::shared_ptr<Subset>& s) {
                return s->potential > newSubset->potential;
            });
            activeSubsets.insert(it, newSubset);
            rootSubsets.remove(s1);
            rootSubsets.remove(s2);
            rootSubsets.push_back(newSubset);

            for (int u : newSubset->vertices) {
                for (int v : newSubset->vertices) {
                    if (u < v) { 
                        auto edge = edges[u][v];

                        const double epsilon = 1e-10;
                        if (active_edges.count(edge) && std::abs(edge->potential) < epsilon) {
                            newSubset->tied_edges.insert(edge);
                        }

                        bool is_erased = active_edges.erase(edge);
                        // std::cout << "active_edges size: " << active_edges.size() << std::endl;
                    }
                }
            }


            // handle tied edges
            if((newSubset->tied_edges).size() > 0) {
                tiedSubsets.push_back(newSubset);
                // cmpTiedEdges(forest, newSubset);
            }

            for(auto s : {s1, s2}) {
                if(!s->is_active) {
                    for(auto edge : active_edges) {
                        if((s->vertices.count(edge->u) + s->vertices.count(edge->v))==1){
                            edge->active_ends += 1;
                        }
                    }
                }
            }
         

        } else {
            // std::cout << "subset go neutral: " << next_subset->potential << std::endl;
            delta_potential = next_subset->potential;
            for (auto& subset : activeSubsets) {
                subset->potential -= delta_potential;
            }
            // for (auto& edge : active_edges) {
            //     edge->potential -= delta_potential * edge->active_ends;
            // }
            updateAndReorderEdges(active_edges, delta_potential);

            // handle subset
            next_subset->is_active = false;

            // handle edge
            for(auto edge : active_edges) {
                if((next_subset->vertices.count(edge->u) + next_subset->vertices.count(edge->v))==1){
                    edge->active_ends -= 1;
                }
            }
            activeSubsets.remove(next_subset);
        }

    }


    return rootSubsets;
}


void prune(std::list<std::shared_ptr<Edge>>& forest, std::shared_ptr<Subset> subset) {
    if (!subset) {
        std::cerr << "Error: subset is null" << std::endl;
        return;
    }

    if (!subset->predecessors.first && !subset->predecessors.second) {
        return;
    }

    if (!subset->is_active) {
        int edgeCount = 0;
        std::list<std::shared_ptr<Edge>> edgesToRemove;

        for (const auto& edge : forest) {
            if (!edge) {
                // std::cerr << "Error: edge is null in forest" << std::endl;
                continue;
            }

            bool delta_edge = (subset->vertices.count(edge->u) + subset->vertices.count(edge->v)) == 1;
            if (delta_edge) {
                edgeCount++;
                edgesToRemove.push_back(edge);
            }
        }

        if (edgeCount == 1) {
            for (const auto& edge : edgesToRemove) {
                forest.remove(edge);
            }
        }
    }

    prune(forest, subset->predecessors.first);
    prune(forest, subset->predecessors.second);
}


void pruneAll(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets) {
    for (auto& s : subsets) {
        if (s) { 
            prune(forest, s);
        }
    }
}

std::unordered_set<int> findLargestConnectedComponent(const std::list<std::shared_ptr<Edge>>& forest) {
    std::unordered_set<int> visited;
    std::unordered_set<int> largestComponent;
    
    for (const auto& edge : forest) {
        if (!edge) continue;

        std::unordered_set<int> currentComponent;
        std::queue<int> queue;

        if (visited.count(edge->u) == 0) {
            queue.push(edge->u);
            visited.insert(edge->u);

            while (!queue.empty()) {
                int node = queue.front();
                queue.pop();
                currentComponent.insert(node);

                for (const auto& e : forest) {
                    if (!e) continue;

                    if (e->u == node && visited.count(e->v) == 0) {
                        queue.push(e->v);
                        visited.insert(e->v);
                    } else if (e->v == node && visited.count(e->u) == 0) {
                        queue.push(e->u);
                        visited.insert(e->u);
                    }
                }
            }

            if (currentComponent.size() > largestComponent.size()) {
                largestComponent = currentComponent;
            }
        }
    }

    return largestComponent;
}


//remark: backtracking, don't parallelize!
int revisitTiedEdges(std::list<std::shared_ptr<Edge>>& forest, std::vector<std::shared_ptr<Subset>>& tiedSubsets, 
                        std::list<std::shared_ptr<Subset>>& rootSubsets, size_t i) {

    if(i == tiedSubsets.size()) {
            std::list<std::shared_ptr<Edge>> forest_(forest.begin(), forest.end());
            pruneAll(forest, rootSubsets);
            std::unordered_set<int> cp = findLargestConnectedComponent(forest);
            return cp.size();
    }

    auto subset = tiedSubsets[i];
    forest.remove(subset->merge_edge);
    int largest_kernel = 0;
    std::shared_ptr<Edge> chosen_edge = nullptr;
    for(auto edge : tiedSubsets[i]->tied_edges) {
        forest.push_back(edge);
        int cp = revisitTiedEdges(forest, tiedSubsets, rootSubsets, i+1);
        if(cp >= largest_kernel) {
            largest_kernel = cp;
            chosen_edge = edge;
        }
        forest.remove(edge);
    }
    forest.push_back(chosen_edge);
    subset->merge_edge = chosen_edge;    
    return largest_kernel;
}



// void revisitTiedEdges(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& tiedSubsets, 
//                         std::list<std::shared_ptr<Subset>>& rootSubsets) {

//     for(auto subset : tiedSubsets) {
//         if (!subset) 
//             continue;
//         if(subset->tied_edges.size() == 0)
//             continue;
//         forest.remove(subset->merge_edge);
//         std::list<std::shared_ptr<Edge>> forest_(forest.begin(), forest.end());

//         subset->tied_edges.insert(subset->merge_edge);
//         int largest_kernel = 0;
//         std::shared_ptr<Edge> chosen_edge = nullptr;
//         for(auto edge : subset->tied_edges) {
//             forest_.push_back(edge);
//             pruneAll(forest_, rootSubsets);
//             std::unordered_set<int> cp = findLargestConnectedComponent(forest_);
//             int kernel_size = cp.size();
//             // std::cout << "Edge: " << edge->u << "-" << edge->v << std::endl;
//             // std::cout << "Largest Kernel: " << kernel_size << std::endl;
//             if(kernel_size > largest_kernel) {
//                 largest_kernel = kernel_size;
//                 chosen_edge = edge;
//             }
//             forest_.remove(edge);
//         }
//         forest.push_back(chosen_edge);
//         subset->merge_edge = chosen_edge;    
//     }

// }

std::shared_ptr<Subset> findSmallestSubsetContainingComponent(std::shared_ptr<Subset> s, 
                                        const std::unordered_set<int>& largestComponent) {
    if (!s) return nullptr;

    bool containsAll = true;
    for (int v : largestComponent) {
        if (s->vertices.count(v) == 0) {
            containsAll = false;
            break;
        }
    }

    if (containsAll) {
        auto leftSubset = findSmallestSubsetContainingComponent(s->predecessors.first, largestComponent);
        auto rightSubset = findSmallestSubsetContainingComponent(s->predecessors.second, largestComponent);

        if (leftSubset) return leftSubset;
        if (rightSubset) return rightSubset;
        return s;
    }

    return nullptr;
}



// void cmpTiedEdges(std::list<std::shared_ptr<Edge>>& forest, std::shared_ptr<Subset> subset) {
//     std::cout << "breaking tie..." << std::endl;
//     if (!subset) {
//         std::cerr << "Error: subset is nullptr" << std::endl;
//         return;
//     }

//     if(subset->tied_edges.size() == 0) {
//         return;
//     }
//     forest.remove(subset->merge_edge);
//     subset->tied_edges.insert(subset->merge_edge);

//     int largest_kernel = 0;
//     std::shared_ptr<Edge> chosen_edge = nullptr;
//     for(auto edge : subset->tied_edges) {
//         if (!edge) {
//             std::cerr << "Error: tied_edges contains a null pointer" << std::endl;
//             continue;
//         }
//         std::list<std::shared_ptr<Edge>> forest_(forest.begin(), forest.end());
//         forest_.push_back(edge);
//         prune(forest_, subset);
//         int kernel_size = LargestKernel(forest_, subset);
//         if(kernel_size > largest_kernel) {
//             largest_kernel = kernel_size;
//             chosen_edge = edge;
//         }
//     }
//     forest.push_back(chosen_edge);
//     subset->merge_edge = chosen_edge;
// }


int LargestKernel(std::list<std::shared_ptr<Edge>>& forest, std::shared_ptr<Subset> s) {
    if (!s) {
        std::cerr << "Error: Subset pointer is null in LargestKernel" << std::endl;
        return 0;
    }
    if (s->vertices.empty()) {
        std::cerr << "Error: Subset has no vertices in LargestKernel" << std::endl;
        return 0;
    }

    std::vector<int> vertices(s->vertices.begin(), s->vertices.end());
    int numVertices = vertices.size();
    UnionFind uf(vertices);

    for (const auto& edge : forest) {
        if (!edge) {
            // std::cerr << "Error: Edge is null in LargestKernel" << std::endl;
            continue;
        }
        if(s->vertices.count(edge->u) + s->vertices.count(edge->v) != 2)
            continue;
        if (uf.find(edge->u) != uf.find(edge->v)) {
            uf.unite(edge->u, edge->v);
        }
    }

    int largestKernelSize = 0;
    for (int i = 0; i < numVertices; ++i) {
        largestKernelSize = std::max(largestKernelSize, uf.getSize(i));
    }
    return largestKernelSize;
}

int findLargestKernel(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets) {
    int largest_kernel = 0;
    std::shared_ptr<Subset> largest_subset = nullptr;

    for(auto s : subsets) {
        int kernel_size = LargestKernel(forest, s);
        if(kernel_size > largest_kernel){
            largest_kernel = kernel_size;
            largest_subset = s;
        }
    }
    return largest_kernel;
}


// std::unordered_set<int> getLargestUnion(UnionFind& uf, const std::shared_ptr<Subset> s) {
//     if (!s) {
//         std::cerr << "Error: subset is null in getLargestUnion" << std::endl;
//         return {};
//     }

//     if (s->vertices.empty()) {
//         std::cerr << "Warning: subset->vertices is empty" << std::endl;
//         return {};
//     }

//     int largestRoot = -1;
//     int maxSize = 0;

//     for (int v : s->vertices) {
//         if (!uf.parent.count(v)) {
//             std::cerr << "Error: vertex " << v << " not found in UnionFind" << std::endl;
//             continue;
//         }

//         int root = uf.find(v);
//         int size = uf.getSize(root);
//         if (size > maxSize) {
//             maxSize = size;
//             largestRoot = root;
//         }
//     }

//     std::unordered_set<int> largestUnion;
//     for (int v : s->vertices) {
//         if (uf.find(v) == largestRoot) {
//             largestUnion.insert(v);
//         }
//     }

//     return largestUnion;
// }


std::vector<std::shared_ptr<Edge>> buildMST(std::list<std::shared_ptr<Edge>>& forest, std::vector<int>& vertices) {
    if (vertices.empty()) {
        std::cerr << "Warning: vertices is empty in buildMST" << std::endl;
        return {};
    }
    std::vector<std::shared_ptr<Edge>> mst;
    std::unordered_set<int> treeVertices(vertices.begin(), vertices.end());

    for (const auto& edge : forest) {
        if (!edge) continue;

        if (treeVertices.count(edge->u) + treeVertices.count(edge->v) == 2) {
            mst.push_back(edge);
        }
    }

    return mst;
}

bool hasIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2) {
    const auto& smaller = (s1.size() < s2.size()) ? s1 : s2;
    const auto& bigger  = (s1.size() < s2.size()) ? s2 : s1;

    for (int x : smaller) {
        if (bigger.find(x) != bigger.end()) {
            return true;
        }
    }
    return false;
}


std::shared_ptr<Subset> pruneSunset(std::shared_ptr<Subset> s, const std::unordered_set<int>& largestComponent) {
    if(!s)
        return s;

    if(hasIntersection(largestComponent, s->vertices)) {
        s->predecessors.first = pruneSunset(s->predecessors.first, largestComponent);
        s->predecessors.second = pruneSunset(s->predecessors.second, largestComponent);
    
        if(s->predecessors.first && s->predecessors.second) {
            s->vertices.clear();
            s->vertices.insert(s->predecessors.first->vertices.begin(), s->predecessors.first->vertices.end());
            s->vertices.insert(s->predecessors.second->vertices.begin(), s->predecessors.second->vertices.end());
            return s;
        }
        else if(s->predecessors.first && !s->predecessors.second)
            return s->predecessors.first;
        else if(!s->predecessors.first && s->predecessors.second)
            return s->predecessors.second;
        else
            return s;
    }
    else {
        return nullptr;
    }

}


void pick(std::shared_ptr<Subset> s, int v, int k, std::vector<int>& pickedVertices) {
    if (!s || k==0) {
        return;
    }

    if(!s->vertices.count(v)) {
        std::cerr << "vertix: " << v << " not in the subset!";
        return;
    }

    auto s1 = s->predecessors.first;
    auto s2 = s->predecessors.second;

    if (!s1 && !s2) {
        pickedVertices.push_back(v);
        return;
    }

    int v1 = s->merge_edge->u;
    int v2 = s->merge_edge->v;

    if (!s1 && s2) {
        pick(s2, v, k, pickedVertices);
        return;
    } else if (!s2 && s1) {
        pick(s1, v, k, pickedVertices);
        return;
    }

    if (!s1->vertices.count(v)) {
        std::swap(s1, s2);
        std::swap(v1, v2);
    }

    int s1_size = s1->vertices.size();

    if (s1_size > k) {
        pick(s1, v, k, pickedVertices);
    } 
    else {
        pickedVertices.insert(pickedVertices.end(), s1->vertices.begin(), s1->vertices.end());
        int remainingK = k - s1_size;
        pick(s2, v2, remainingK, pickedVertices);
    } 
}


std::vector<std::shared_ptr<Edge>> FindMST(std::list<std::shared_ptr<Edge>>& forest, std::list<std::shared_ptr<Subset>>& subsets, int k) {
    int largest_kernel = 0;
    std::shared_ptr<Subset> largest_subset = nullptr;

    for (const auto& s : subsets) {
        if (!s) {
            std::cerr << "Warning: null subset in FindMST" << std::endl;
            continue;
        }

        int kernel_size = LargestKernel(forest, s);
        if (kernel_size > largest_kernel) {
            largest_kernel = kernel_size;
            largest_subset = s;
        }
    }

    if (!largest_subset) {
        std::cerr << "Error: no valid largest_subset found in FindMST" << std::endl;
        return {};
    }

    std::vector<int> pickedVertices;
    pick(largest_subset, *largest_subset->vertices.begin(), k, pickedVertices);
    std::vector<std::shared_ptr<Edge>>mst = buildMST(forest, pickedVertices);
    return mst;

}


double findLambda1(const std::vector<std::vector<double>>& costs, int k, double epsilon) {
    double lambda_min = 0;
    double lambda_max = 0;

    for (size_t i = 0; i < costs.size(); ++i) {        
        for (size_t j = i+1; j < costs.size(); ++j) {       
            if (costs[i][j] > lambda_max) {     
                lambda_max = costs[i][j];
            }
        }
    }

    while (lambda_max - lambda_min > epsilon) {
        // std::cout << "lambda_max: " << lambda_max << std::endl;
        // std::cout << "lambda_min: " << lambda_min << std::endl;
        double lambda = lambda_min + (lambda_max - lambda_min) / 2.0;

        std::list<std::shared_ptr<Edge>> forest;
        std::list<std::shared_ptr<Subset>> tiedSubsets;
        auto subsets = primalDualSubroutine(costs, lambda, forest, tiedSubsets);
        std::vector<std::shared_ptr<Subset>> tiedSubsets_(tiedSubsets.begin(), tiedSubsets.end());
        revisitTiedEdges(forest, tiedSubsets_, subsets, 0);
        pruneAll(forest, subsets);
        int largest_kernel = findLargestKernel(forest, subsets);

        if (largest_kernel >= k) {
            lambda_max = lambda;
        } else {
            lambda_min = lambda;
        }
    }
    return lambda_max;
}


std::vector<std::shared_ptr<Edge>> kMST(const std::vector<std::vector<double>>& costs, int k) {

    std::list<std::shared_ptr<Edge>> forest;
    std::list<std::shared_ptr<Subset>> tiedSubsets;

    double lambda = findLambda1(costs, k, 1e-9);
    // std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Find Lambda: " << lambda << std::endl;

    std::list<std::shared_ptr<Subset>> subsets = primalDualSubroutine(costs, lambda, forest, tiedSubsets);
    // revisitTiedEdges(forest, tiedSubsets, subsets);
    std::vector<std::shared_ptr<Subset>> tiedSubsets_(tiedSubsets.begin(), tiedSubsets.end());
    revisitTiedEdges(forest, tiedSubsets_, subsets, 0);
    pruneAll(forest, subsets);
    // std::cout << "num root subsets: " << subsets.size() << std::endl;
    // for(auto s:subsets) {
    //     std::cout << s->vertices.size() << ", ";
    // }
    // std::cout  << std::endl;

    std::unordered_set<int> largestComponent = findLargestConnectedComponent(forest);
    // std::cout << "largestComponent: " << largestComponent.size() << std::endl;
    std::shared_ptr<Subset> smallestSubset = nullptr;

    for(auto s: subsets) {
        smallestSubset = findSmallestSubsetContainingComponent(s, largestComponent);
        if (smallestSubset)
            break;
    }
    // std::cout << "smallestSubset: " << smallestSubset->vertices.size() << std::endl;

    std::unordered_set<std::shared_ptr<Edge>> Forest(forest.begin(), forest.end());
    smallestSubset = pruneSunset(smallestSubset, largestComponent);
    // if (!smallestSubset) {
    //     std::cout << "smallestSubset is null after prune!" << std::endl;
    // } else {
    //     std::cout << "smallestSubset after prune: " 
    //             << smallestSubset->vertices.size() << std::endl;
    // }

    // std::cout << "smallestSubset after prune: " << smallestSubset->vertices.size() << std::endl;


    std::vector<int> pickedVertices;
    pick(smallestSubset, *smallestSubset->vertices.begin(), k, pickedVertices);
    // std::cout << "pickedVertices: " << pickedVertices.size() << std::endl;
    std::vector<std::shared_ptr<Edge>>mst = buildMST(forest, pickedVertices);

    // std::vector<std::shared_ptr<Edge>> mst = FindMST(forest, subsets, k);
    // double total_weight = 0;
    // std::cout << "tree size: " << mst.size() << std::endl;
    // for(auto e: mst) {
    //     // std::cout << e->u << "-" << e->v << ": " << e->weight << ", ";
    //     total_weight += e->weight;
    // }
    // std::cout << std::endl << "total_weight: " << total_weight << std::endl;

    return mst;

}


std::vector<int> findK(const std::vector<std::vector<double>>& costs, double budget) {
    std::vector<int> best_path;
    double best_cost = std::numeric_limits<double>::infinity();

    int k_min = 1, k_max = costs.size();
    while (k_min <= k_max) {
        int k = k_min + (k_max - k_min) / 2;

        std::vector<std::shared_ptr<Edge>> mst = kMST(costs, k);
        std::vector<int> path = christofides(costs, mst, mst[0]->u);

        double path_cost = 0.0;
        for (size_t i = 0; i < path.size() - 1; i++) {
            path_cost += costs[path[i]][path[i + 1]];
        }

        std::cout << "k: " << k << ", Path Cost: " << path_cost << std::endl;

        if (path_cost <= budget) {
            if (path.size() > best_path.size() || (path.size() == best_path.size() && path_cost < best_cost)) {
                best_cost = path_cost;
                best_path = path;
            }
            k_min = k + 1;
        } else {
            k_max = k - 1;
        }
    }

    return best_path;
}

void findBestTour(const std::vector<std::vector<double>>& costs, double budget) {
    std::vector<int>best_tour =  findK(costs, budget);
    std::cout << " Best Prize: " << best_tour.size() << std::endl;
    std::cout << "Best Tour: ";
    for (int city : best_tour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    double total_cost = 0.0;
    for (size_t i = 0; i < best_tour.size() - 1; i++) {
        total_cost += costs[best_tour[i]][best_tour[i + 1]];
        // std::cout << total_cost << " ";
    }
    std::cout << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;
}

