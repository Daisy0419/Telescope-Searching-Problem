#include "ReadData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>

// // in degree
// double angularDistance(double ra1, double dec1, double ra2, double dec2) {
//     return acos(sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1 - ra2)) * (180.0 / M_PI);
// }

#define DEG_TO_RAD (M_PI / 180.0)

double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG_TO_RAD;
    double dec1_rad = dec1 * DEG_TO_RAD;
    double ra2_rad = ra2 * DEG_TO_RAD;
    double dec2_rad = dec2 * DEG_TO_RAD;

    // Difference within (-π, π]
    double delta_ra = fmod(ra1_rad - ra2_rad + M_PI, 2 * M_PI) - M_PI;

    double angle_rad = acos(sin(dec1_rad) * sin(dec2_rad) +
                            cos(dec1_rad) * cos(dec2_rad) * cos(delta_ra));

    return angle_rad * (180.0 / M_PI);
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
    getline(file, line);  // Skip header

    costs.clear();
    probability.clear();  // Ensure no leftover values

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
    costs.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                costs[i][j] = angularDistance(ra[i], dec[i], ra[j], dec[j]) / telescope_speed + dwell_time;
                costs[j][i] = costs[i][j];
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
}


void read_data_filtered(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<double>& ranks, 
               double telescope_speed, double dwell_time, double select_probability) {
    std::vector<double> ra, dec;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    costs.clear();
    probability.clear();
    costs_filtered.clear();
    probability_filtered.clear();
    ranks.clear();
    ra.clear();
    dec.clear();

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ranks.push_back(std::stod(row[0])); 
        ra.push_back(std::stod(row[2]));  
        dec.push_back(std::stod(row[3]));  
        probability.push_back(std::stod(row[4]));  
    }
    file.close();

    // Create a vector of indices for sorting
    size_t n = probability.size();
    std::vector<size_t> indices(n);
    for (size_t i = 0; i < n; ++i) {
        indices[i] = i;
    }

    // Sort indices
    std::sort(indices.begin(), indices.end(), 
              [&probability](size_t a, size_t b) { return probability[a] > probability[b]; });

    // Reorder all vectors according to the sorted indices
    std::vector<double> temp_probability(n), temp_ranks(n), temp_ra(n), temp_dec(n);
    for (size_t i = 0; i < n; ++i) {
        temp_probability[i] = probability[indices[i]];
        temp_ranks[i] = ranks[indices[i]];
        temp_ra[i] = ra[indices[i]];
        temp_dec[i] = dec[indices[i]];
    }

    probability = std::move(temp_probability);
    ranks = std::move(temp_ranks);
    ra = std::move(temp_ra);
    dec = std::move(temp_dec);

    double sum_probability = 0.0;
    int select_num = n;
    for (size_t i = 0; i < n; ++i) {
        probability_filtered.push_back(probability[i]);
        sum_probability += probability[i];
        if (sum_probability > select_probability) {
            select_num = i + 1;
            break;
        }
    }

    costs_filtered.assign(select_num, std::vector<double>(select_num, 0.0));
    for (size_t i = 0; i < select_num; ++i) {
        for (size_t j = i+1; j < select_num; ++j) {
            if (i != j) {
                costs_filtered[i][j] = angularDistance(ra[i], dec[i], ra[j], dec[j]) / telescope_speed + dwell_time;
                costs_filtered[j][i] = costs_filtered[i][j];  // Symmetry in cost
            }
        }
        costs_filtered[i][i] = std::numeric_limits<double>::infinity(); 
    }

    // size_t n = ra.size();
    costs.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                costs[i][j] = angularDistance(ra[i], dec[i], ra[j], dec[j]) / telescope_speed + dwell_time;
                costs[j][i] = costs[i][j];
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
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

    costs.clear();
    probability.clear();

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
    costs.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                // double travel_time = std::max(std::abs(ra[j] - ra[i]), std::abs(dec[j] - dec[i]));
                double ra_diff = std::abs(ra[j] - ra[i]);
                ra_diff = std::min(ra_diff, 360.0 - ra_diff);
                double dec_diff = std::abs(dec[j] - dec[i]);
                double travel_time = std::max(ra_diff, dec_diff);
                costs[i][j] = (telescope_speed > 0) ? (travel_time / telescope_speed + dwell_time) : std::numeric_limits<double>::infinity();
                costs[j][i] = costs[i][j];
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
}



void read_data_deep_slow_filtered(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<double>& ranks, 
               double telescope_speed, double dwell_time, double select_probability) {
    std::vector<double> ra, dec;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    costs.clear();
    probability.clear();
    costs_filtered.clear();
    probability_filtered.clear();
    ranks.clear();
    ra.clear();
    dec.clear();

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ranks.push_back(std::stod(row[0])); 
        ra.push_back(std::stod(row[2]));  
        dec.push_back(std::stod(row[3]));  
        probability.push_back(std::stod(row[4]));  
    }
    file.close();

    // Create a vector of indices for sorting
    size_t n = probability.size();
    std::vector<size_t> indices(n);
    for (size_t i = 0; i < n; ++i) {
        indices[i] = i;
    }

    // Sort indices
    std::sort(indices.begin(), indices.end(), 
              [&probability](size_t a, size_t b) { return probability[a] > probability[b]; });

    // Reorder all vectors according to the sorted indices
    std::vector<double> temp_probability(n), temp_ranks(n), temp_ra(n), temp_dec(n);
    for (size_t i = 0; i < n; ++i) {
        temp_probability[i] = probability[indices[i]];
        temp_ranks[i] = ranks[indices[i]];
        temp_ra[i] = ra[indices[i]];
        temp_dec[i] = dec[indices[i]];
    }

    probability = std::move(temp_probability);
    ranks = std::move(temp_ranks);
    ra = std::move(temp_ra);
    dec = std::move(temp_dec);

    double sum_probability = 0.0;
    int select_num = n;
    for (size_t i = 0; i < n; ++i) {
        probability_filtered.push_back(probability[i]);
        sum_probability += probability[i];
        if (sum_probability > select_probability) {
            select_num = i + 1;
            break;
        }
    }

    costs_filtered.assign(select_num, std::vector<double>(select_num, 0.0));
    for (size_t i = 0; i < select_num; ++i) {
        for (size_t j = 0; j < select_num; ++j) {
            if (i != j) {
                // double travel_time = std::max(std::abs(ra[j] - ra[i]), std::abs(dec[j] - dec[i]));
                double ra_diff = std::abs(ra[j] - ra[i]);
                ra_diff = std::min(ra_diff, 360.0 - ra_diff);
                double dec_diff = std::abs(dec[j] - dec[i]);
                double travel_time = std::max(ra_diff, dec_diff);
                costs_filtered[i][j] = (telescope_speed > 0) ? (travel_time / telescope_speed + dwell_time) : std::numeric_limits<double>::infinity();
                costs_filtered[j][i] = costs_filtered[i][j];
            }
        }
        costs_filtered[i][i] = std::numeric_limits<double>::infinity(); 
    }

    costs.assign(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                // double travel_time = std::max(std::abs(ra[j] - ra[i]), std::abs(dec[j] - dec[i]));
                double ra_diff = std::abs(ra[j] - ra[i]);
                ra_diff = std::min(ra_diff, 360.0 - ra_diff);
                double dec_diff = std::abs(dec[j] - dec[i]);
                double travel_time = std::max(ra_diff, dec_diff);
                costs[i][j] = (telescope_speed > 0) ? (travel_time / telescope_speed + dwell_time) : std::numeric_limits<double>::infinity();
                costs[j][i] = costs[i][j];
            }
        }
        costs[i][i] = std::numeric_limits<double>::infinity(); 
    }
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