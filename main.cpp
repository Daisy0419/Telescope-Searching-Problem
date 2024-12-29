#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "helplers.h"
#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "Christofides.h"
#include "kMST.h"

std::vector<std::vector<double>> initialize_distances(int NUM_CITIES) {
    std::vector<std::vector<double>> costs(NUM_CITIES, std::vector<double>(NUM_CITIES, std::numeric_limits<double>::infinity()));
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = i + 1; j < NUM_CITIES; j++) {
            costs[i][j] = random_double(10.0, 100.0); 
            costs[j][i] =  costs[i][j];
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
    return costs;
}

std::vector<std::vector<double>> initialize_distances(const std::string &fileName) {
    std::vector<std::vector<double>> costs;
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return costs;
    }

    std::string line, ignoredValue, graphName, type;
    int dimension = 0;

    std::getline(infile, line); // NAME
    std::istringstream issName(line);
    issName >> ignoredValue >> graphName;

    infile >> ignoredValue >> type;        // TYPE
    infile >> ignoredValue >> dimension;  // DIMENSION


    costs.resize(dimension, std::vector<double>(dimension, std::numeric_limits<double>::infinity()));
    // Skip to EDGE_WEIGHT_SECTION
    while (std::getline(infile, line)) {
        if (line == "EDGE_WEIGHT_SECTION") {
            break;
        }
    }

    // Parse edges and populate the adjacency matrix
    int head, tail;
    double distance;
    while (std::getline(infile, line)) {
        std::istringstream edgeStream(line);
        if (!(edgeStream >> head >> tail >> distance)) {
            std::cerr << "Error: Malformed edge line '" << line << "'" << std::endl;
            continue;
        }
        costs[head][tail] = distance;
        costs[tail][head] = distance;
    }
    std::cout << "Distances successfully loaded from " << fileName << std::endl;
    return costs;
}

int main() {
    
    std::string file = "../bier127.tsp";
    double budget = 47358.8;
    std::vector<std::vector<double>> costs = initialize_distances(file);

    // double budget = 500;
    // std::vector<std::vector<double>> costs = initialize_distances(200);

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    // ant colony
    std::cout << "*********running ant colony*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AntColony ac = AntColony(costs);
    ac.ant_colony_optimization(budget);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;



    // genetic
    std::cout << "*********running genetic*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    genetic_optimization(costs, budget, 0);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;


    //greedy
    std::cout << "*********running greedy*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    max_rooted_path(costs, budget, 0);

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    //kmst
    std::cout << "*********running kmst*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    findBestTour(costs, budget);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    return 0;
}

