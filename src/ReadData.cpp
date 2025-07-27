#include "ReadData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)

double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG_TO_RAD;
    double dec1_rad = dec1 * DEG_TO_RAD;
    double ra2_rad = ra2 * DEG_TO_RAD;
    double dec2_rad = dec2 * DEG_TO_RAD;

    double delta_ra = fmod(ra1_rad - ra2_rad + M_PI, 2 * M_PI) - M_PI;

    double cos_angle = sin(dec1_rad) * sin(dec2_rad) +
                       cos(dec1_rad) * cos(dec2_rad) * cos(delta_ra);

    // Clamp to [-1, 1]
    cos_angle = std::min(1.0, std::max(-1.0, cos_angle));

    double angle_rad = acos(cos_angle);

    return angle_rad * RAD_TO_DEG;
}

// read rank, ra, dec and prize form tiling file
void read_data_from_file(const std::string& filename, std::vector<double>& probability, 
                        std::vector<int>& ranks, std::vector<double>& ra, std::vector<double>& dec) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;

        while (getline(ss, token, ',')) {
            row.push_back(token);
        }

        ranks.push_back(std::stoi(row[0])); 
        ra.push_back(std::stod(row[2]));  
        dec.push_back(std::stod(row[3]));  
        probability.push_back(std::stod(row[4]));  
    }
    file.close();

}

// edge cost = geodesic distancing + dwell time (for unified dwelltime)
// cost[i][j] = anguler_distance[i][j]/slew_rate + dwelltime[j]
void compute_costs_geodesic(const std::vector<double>& ra, const std::vector<double>& dec, 
                        std::vector<std::vector<double>>& costs, 
                        double slew_rate, double dwell_time) {

    size_t n = ra.size();
    costs.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                costs[i][j] = angularDistance(ra[i], dec[i], ra[j], dec[j]) / slew_rate + dwell_time;
                costs[j][i] = costs[i][j];
            }
        }
        // costs[i][i] = std::numeric_limits<double>::infinity(); 
        costs[i][i] = dwell_time;
    }
}

// edge cost = max axes distancing + dwelltime
// cost[i][j] = axes_distancing[i][j]/slew_rate + dwelltime[j]
void compute_costs_deepslow(const std::vector<double>& ra, const std::vector<double>& dec, 
                        std::vector<std::vector<double>>& costs, 
                        double slew_rate, double dwell_time) {

    size_t n = ra.size();
    costs.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double ra_diff = std::abs(ra[j] - ra[i]);
                ra_diff = std::min(ra_diff, 360.0 - ra_diff);
                double dec_diff = std::abs(dec[j] - dec[i]);
                double travel_dist = std::max(ra_diff, dec_diff);
                costs[i][j] = (slew_rate > 0) ? (travel_dist / slew_rate + dwell_time) : std::numeric_limits<double>::infinity();
                costs[j][i] = costs[i][j];
            }
        }
        // costs[i][i] = std::numeric_limits<double>::infinity(); 
        costs[i][i] = dwell_time;
    }
}


// Computes air mass using the Kasten & Young model
double computeAirMass(double zenithAngleDegrees) {
    // if (zenithAngleDegrees >= 90.0) {
    //     std::cerr << "Zenith angle must be less than 90 degrees.";
    // }
    double theta = DEG_TO_RAD * zenithAngleDegrees;
    double numerator = 1.0;
    double denominator = std::cos(theta) + 0.50572 * std::pow(96.07995 - zenithAngleDegrees, -1.6364);
    return numerator / denominator;
}

// Computes the light attenuation scaling factor for a given air mass
double computeScalingFactor(double airMass) {
    return 1.1129 * std::exp(-0.107 * airMass);
}

// Computes the scaled dwell time needed to achieve same SNR as reference dwell time t_ref at s=1
double computeDwellTime(double t_ref, double scalingFactor) {
    return t_ref / (scalingFactor * scalingFactor);
}


// Orienteering Graph
// cost[i][j] = anguler_distance[i][j]/slew_rate + 0.5*dwelltime[i] + 0.5*dwelltime[j]
void compute_costs_orienteering(const std::vector<double>& ra, const std::vector<double>& dec, 
                        const std::vector<double>& dwell_times,
                        std::vector<std::vector<double>>& costs, 
                        double slew_rate, int end_pos, double& _padding) {

    size_t n = ra.size();
    if(n != end_pos) {
        std::cerr << "Wrong sizing in compute_costs_orienteering func" << std::endl;
    }

    costs.assign(n + 1, std::vector<double>(n + 1, 0.0)); //to include end_pos t'

    _padding = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (i != j) {
                double slew_time = angularDistance(ra[i], dec[i], ra[j], dec[j]) / slew_rate;
                costs[i][j] = slew_time + 0.5 * dwell_times[i] + 0.5 * dwell_times[j];
                costs[j][i] = costs[i][j];
                if(slew_time > _padding) _padding = slew_time;
            }
        }
        // costs[i][i] = std::numeric_limits<double>::infinity(); 
        costs[i][i] = dwell_times[i];
    }
    _padding = ceil(_padding); //smallest big enough padding
    for (size_t i = 0; i < n; ++i) {
            costs[i][end_pos] = 0.5 * dwell_times[i] + _padding;
            costs[end_pos][i] = costs[i][end_pos];
    }

}


// build graph for s-t orienteering
std::tuple<int, int, double> buildGraphOrienteering(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks, std::vector<double>& dwell_times,
               double slew_rate, bool is_deepslow, int init_pos) {
    std::vector<double> ra, dec;

    // get probability, ranks, ra, dec
    read_data_from_file(filename, probability, ranks, ra, dec);
    double init_ra, init_dec;
    if (init_pos >= 0 && init_pos < (int)ra.size()) {
        init_ra = ra[init_pos], init_dec = dec[init_pos];
    } else {
        throw std::invalid_argument("Invalid init position.");
    }
    int end_rank = probability.size();

    dwell_times.resize(probability.size());
    double dwelltime_ref = 1.0;
    double ra_zenith = ra[0], dec_zenith = dec[0];

    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < probability.size(); ++i) {
        double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
        if (zenithAngle >= 90.0) {
            continue; // Skip this index
        }
        double airMass = computeAirMass(zenithAngle);
        double scalingFactor = computeScalingFactor(airMass);
        dwell_times[i] = computeDwellTime(dwelltime_ref, scalingFactor);
        valid_indices.push_back(i);
    }

    // Filter all vectors based on valid_indices
    std::vector<double> filtered_ra, filtered_dec, filtered_probability, filtered_dwell_times;
    std::vector<int> filtered_ranks;
    for (size_t idx : valid_indices) {
        filtered_ra.push_back(ra[idx]);
        filtered_dec.push_back(dec[idx]);
        filtered_probability.push_back(probability[idx]);
        filtered_ranks.push_back(ranks[idx]);
        filtered_dwell_times.push_back(dwell_times[idx]);
    }
    ra = std::move(filtered_ra);
    dec = std::move(filtered_dec);
    probability = std::move(filtered_probability);
    ranks = std::move(filtered_ranks);
    dwell_times = std::move(filtered_dwell_times);

    int start_pos = probability.size();
    ra.push_back(init_ra);
    dec.push_back(init_dec);
    probability.push_back(0.0);
    ranks.push_back(0);
    dwell_times.push_back(0.0);

    // Push as end node
    int end_pos = probability.size();
    probability.push_back(0.0);
    ranks.push_back(end_rank);
    dwell_times.push_back(0.0);

    double padding_value;
    compute_costs_orienteering(ra, dec, dwell_times, costs, slew_rate, end_pos, padding_value);

    return {start_pos, end_pos, padding_value};
}

// rebuild graph for s-t orienteering to exlude a existing path 
// and set the last node in the existing path as the start node
std::pair<int, int> removeNodesOrienteering(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& dwell_times,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping,
                int start_idx, int end_idx) {

    costs_removed.clear();
    prizes_removed.clear();
    ranks_mapping.clear();

    if (nodes_to_remove.empty()) {
        costs_removed = costs;
        prizes_removed = prizes;
        ranks_mapping.resize(prizes.size());
        std::iota(ranks_mapping.begin(), ranks_mapping.end(), 0);
        return {start_idx, end_idx};
    }
        
    std::unordered_set<int> remove_set(nodes_to_remove.begin(), nodes_to_remove.end());
    remove_set.erase(start_idx);
    remove_set.erase(end_idx);

    std::vector<int> kept;
    kept.reserve(prizes.size() - remove_set.size());

    kept.push_back(start_idx);
    kept.push_back(end_idx);

    for (int old_idx = 0; old_idx < int(prizes.size()); ++old_idx)
        if (remove_set.find(old_idx) == remove_set.end() && old_idx != start_idx && old_idx != end_idx)
            kept.push_back(old_idx);


    const std::size_t S = kept.size();

    costs_removed.assign(S, std::vector<double>(S));
    prizes_removed.resize(S);

    ranks_mapping = kept;

    for (std::size_t i = 0; i < S; ++i) {
        const int old_i = kept[i];
        prizes_removed[i] = prizes[old_i];
        for (std::size_t j = i+1; j < S; ++j) {
            const int old_j = kept[j];
            costs_removed[i][j] = costs[old_i][old_j];
            if(old_i == start_idx || old_j == start_idx) 
                costs_removed[i][j] -= 0.5 * dwell_times[start_idx];
            costs_removed[j][i] = costs_removed[i][j];
        }
    }

    int new_start_idx = 0;
    int new_end_idx = 1;
    return {new_start_idx, new_end_idx}; 
}
