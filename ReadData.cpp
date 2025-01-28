#include "ReadData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <cmath>

// in degree
double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    return acos(sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1 - ra2)) * (180.0 / M_PI);
}

void read_data(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, double telescope_speed, double dwell_time) {
    std::vector<double> ra, dec;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line); 

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ra.push_back(std::stod(row[2])); 
        dec.push_back(std::stod(row[3])); 
        probability.push_back(std::stod(row[4]));
    }
    file.close();

    size_t n = ra.size();
    costs.resize(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                costs[i][j] = angularDistance(ra[i], dec[i], ra[j], dec[j]) / telescope_speed + dwell_time;
                costs[j][i] = costs[i][j];
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }

    // std::cout << "Distance Matrix (Angular Distances in Degrees):" << std::endl;
    // for (const auto& row : costs) {
    //     for (double val : row) {
    //         std::cout << std::fixed << std::setprecision(6) << val << "\t";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "\nPrize Vector (Probabilities):" << std::endl;
    // for (double val : probability) {
    //     std::cout << std::fixed << std::setprecision(6) << val << " ";
    // }
    // std::cout << std::endl;
}

void read_data_deep_slow(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, double telescope_speed, double dwell_time) {
    std::vector<double> ra, dec;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line); 

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ra.push_back(std::stod(row[2])); 
        dec.push_back(std::stod(row[3])); 
        probability.push_back(std::stod(row[4]));
    }
    file.close();

    size_t n = ra.size();
    costs.resize(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double travel_time = std::max(std::abs(ra[j] - ra[i]), std::abs(dec[j] - dec[i])) / telescope_speed;;
                costs[i][j] = travel_time + dwell_time;
                costs[j][i] = costs[i][j]; 
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }

}

//build random graph
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

//build graph from tsp data
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