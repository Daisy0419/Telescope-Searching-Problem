
#include <fstream>
#include <sstream>
#include <random>

#include "graph.h"
#include "helplers.h"

//read bier127.tsp
Graph::Graph(const std::string &fileName) {
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
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
}


//random generate
Graph::Graph(int NUM_CITIES) {
    costs.resize(NUM_CITIES, std::vector<double>(NUM_CITIES, std::numeric_limits<double>::infinity()));
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = i + 1; j < NUM_CITIES; j++) {
            costs[i][j] = random_double(10.0, 100.0); 
            costs[j][i] =  costs[i][j];
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
}
