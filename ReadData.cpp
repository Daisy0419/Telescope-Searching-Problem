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

// read rank, ra, dec and prize form file
void read_data_from_file(const std::string& filename, std::vector<double>& probability, 
                        std::vector<int>& ranks, std::vector<double>& ra, std::vector<double>& dec) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line, token;
    getline(file, line);  // Skip header

    // probability.clear();
    // ranks.clear();
    // ra.clear();
    // dec.clear();

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

    costs.assign(n + 1, std::vector<double>(n + 1, 0.0));//to include end_pos t'

    _padding = 0.0;
    // std::cout << "slew_rate: " << slew_rate << std::endl;
    // std::cout << "slew_time: ";
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
    // std::cout << "_padding: " << _padding << std::endl;
    _padding = ceil(_padding); //smallest big enough padding
    // std::cout << "_padding: " << _padding << std::endl;
    for (size_t i = 0; i < n; ++i) {
            costs[i][end_pos] = 0.5 * dwell_times[i] + _padding;
            costs[end_pos][i] = costs[i][end_pos];
    }

}

std::pair<double, double> geometricCenter(const std::vector<double>& ra_deg,
                                               const std::vector<double>& dec_deg) {
    double x = 0.0, y = 0.0, z = 0.0;

    for (size_t i = 0; i < ra_deg.size(); ++i) {
        double ra_rad = ra_deg[i] * DEG_TO_RAD;
        double dec_rad = dec_deg[i] * DEG_TO_RAD;

        x += std::cos(dec_rad) * std::cos(ra_rad);
        y += std::cos(dec_rad) * std::sin(ra_rad);
        z += std::sin(dec_rad);
    }

    x /= ra_deg.size();
    y /= ra_deg.size();
    z /= ra_deg.size();

    double hyp = std::sqrt(x * x + y * y);
    double ra_center = std::atan2(y, x);         // in radians
    double dec_center = std::atan2(z, hyp);      // in radians

    if (ra_center < 0) ra_center += 2 * M_PI;

    return { ra_center * RAD_TO_DEG, dec_center * RAD_TO_DEG };
}

// std::tuple<int, int, double> buildGraphOrienteering(const std::string& filename, std::vector<std::vector<double>>& costs, 
//                std::vector<double>& probability, std::vector<int>& ranks,
//                double slew_rate, bool is_deepslow, int init_pos) {
//     std::vector<double> ra, dec;

//     //get probability, ranks, ra, dec
//     read_data_from_file(filename, probability, ranks, ra, dec);

//     //assume dwell times reverse to probability
//     std::vector<double> dwell_times(probability.size());
//     double dwelltime_ref = 1.0;
//     double ra_zenith = ra[0], dec_zenith = dec[0];
//     // auto [ra_zenith, dec_zenith] = geometricCenter(ra, dec);

//     for (size_t i = 0; i < probability.size(); ++i) {
//         double zenithAngle = angularDistance(ra[i], dec[i], ra_zenith, dec_zenith);
//         if(zenithAngle >= 90.0) {
//             // dwell_times[i] = std::numeric_limits<double>::infinity();
//             dwell_times[i] = 10E6;
//             probability[i] = 0.0;
//             continue;
//         }
//         double airMass = computeAirMass(zenithAngle);
//         double scalingFactor = computeScalingFactor(airMass);
//         dwell_times[i] = computeDwellTime(dwelltime_ref, scalingFactor);
//     }

//     int start_pos = probability.size();
//     if (init_pos >= 0 && init_pos < (int)ra.size()) {
//         double init_ra = ra[init_pos], init_dec = dec[init_pos];
//         ra.push_back(init_ra);
//         dec.push_back(init_dec);
//         probability.push_back(0.0);
//         ranks.push_back(0);
//         dwell_times.push_back(0.0);
//     }

//     // Push as end node
//     int end_pos = probability.size();
//     probability.push_back(0.0);
//     ranks.push_back(end_pos);
//     dwell_times.push_back(0.0);

//     double padding_value;
//     compute_costs_orienteering(ra, dec, dwell_times, costs, slew_rate, end_pos, padding_value);

//     return {start_pos, end_pos, padding_value};
// }


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

    // assume dwell times reverse to probability
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

//init_pos in front
void sortNodes(const std::vector<double>& prizes, std::vector<double>& prizes_sorted, 
                const std::vector<int>& ranks, std::vector<int>& ranks_sorted, 
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
    ranks_sorted.resize(n);
    indices.resize(n);
    
    prizes_sorted[0] = prizes[init_pos];
    ranks_sorted[0] = ranks[init_pos];
    indices[0] = init_pos;

    size_t j = 1;
    for (size_t i = 0; i < n; ++i) {
        if (indices_temp[i] == init_pos)
            continue;
        prizes_sorted[j] = prizes[indices_temp[i]];
        ranks_sorted[j] = ranks[indices_temp[i]];
        indices[j] = indices_temp[i];
        ++j;
    }
}


// filter the regions with sum probability of "select_probability"
void chooseTopRegion(const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes_sorted, 
                    const std::vector<int>& ranks_sorted, 
                    const std::vector<int>& indices,
                    std::vector<std::vector<double>>& costs_chosenN, 
                    std::vector<double>& prizes_chosenN,
                    std::vector<int>& ranks_chosenN,
                    const double target_prize) {

    if (target_prize < 0) {
        costs_chosenN.clear();
        ranks_chosenN.clear();
        prizes_chosenN.clear();
        return;
    }

    double cur_prize = 0.0;
    int N = prizes_sorted.size();
    int ChosenN = 0;
    for (; ChosenN < N; ++ChosenN) {
        cur_prize += prizes_sorted[ChosenN];
        if (cur_prize >= target_prize) {
            ++ChosenN;
            break;
        }
    }

    ranks_chosenN.assign(ranks_sorted.begin(), ranks_sorted.begin() + ChosenN);
    prizes_chosenN.assign(prizes_sorted.begin(), prizes_sorted.begin() + ChosenN);

    costs_chosenN.resize(ChosenN, std::vector<double>(ChosenN));
    for (int i = 0; i < ChosenN; ++i) {
        for (int j = 0; j < ChosenN; ++j) {
            costs_chosenN[i][j] = costs[indices[i]][indices[j]];
        }
    }
}


//return new init_pos, which set to be the last node in nodes_to_remove
int removeNodes(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping) {

    costs_removed.clear();
    prizes_removed.clear();
    ranks_mapping.clear();

    if (nodes_to_remove.empty()) {
        costs_removed = costs;
        prizes_removed = prizes;
        ranks_mapping.resize(prizes.size());
        std::iota(ranks_mapping.begin(), ranks_mapping.end(), 0);
        return 0;
    }
        
    std::unordered_set<int> remove_set(nodes_to_remove.begin(), nodes_to_remove.end());

    std::vector<int> kept;
    kept.reserve(prizes.size() - remove_set.size() + 1);

    const int root_old_idx = nodes_to_remove.back();
    kept.push_back(root_old_idx);

    for (int old_idx = 0; old_idx < static_cast<int>(prizes.size()); ++old_idx)
        if (remove_set.find(old_idx) == remove_set.end() && old_idx != root_old_idx)
            kept.push_back(old_idx);

    const std::size_t S = kept.size();

    costs_removed.assign(S, std::vector<double>(S));
    prizes_removed.resize(S);

    ranks_mapping = kept;

    for (std::size_t i = 0; i < S; ++i)
    {
        const int old_i = kept[i];
        prizes_removed[i] = prizes[old_i];

        for (std::size_t j = 0; j < S; ++j)
        {
            const int old_j = kept[j];
            costs_removed[i][j] = costs[old_i][old_j];
        }
    }
    return 0; 
}


//pass init_pos in form {ra, dec}
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow,
               std::pair<double, double> init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (!std::isnan(init_pos.first) && !std::isnan(init_pos.second)) {
        ra.push_back(init_pos.first);
        dec.push_back(init_pos.second);
        probability.push_back(0.0);
        ranks.push_back(0);
    }

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }

    int start_pos = probability.size() - 1;
    return start_pos;
}

//pass init_pos in form index
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow,
               int init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (init_pos >= 0 && init_pos < (int)ra.size()) {
        double init_ra = ra[init_pos], init_dec = dec[init_pos];
        ra.push_back(init_ra);
        dec.push_back(init_dec);
        probability.push_back(0.0);
        ranks.push_back(0);
    }

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }

    int start_pos = probability.size() - 1;
    return start_pos;
}


//pass init_pos in form {ra, dec}; filtere region with sum prob of target_sum_prob
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               std::pair<double, double> init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (!std::isnan(init_pos.first) && !std::isnan(init_pos.second)) {
        ra.push_back(init_pos.first);
        dec.push_back(init_pos.second);
        probability.push_back(0.0);
        ranks.push_back(0);
    }

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }

    std::vector<double> probability_sorted;
    std::vector<int> ranks_sorted;
    std::vector<int> indices;
    int start_idx = probability.size() - 1;
    sortNodes(probability, probability_sorted, ranks, ranks_sorted, indices, start_idx);
    chooseTopRegion(costs, probability_sorted, ranks_sorted, indices, costs_filtered, probability_filtered, ranks_filtered, target_sum_prob);

    start_idx = 0;  // after choosing top region, initial point is at index 0
    return start_idx;
}

int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               int init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (init_pos >= 0 && init_pos < (int)ra.size()) {
        double init_ra = ra[init_pos], init_dec = dec[init_pos];
        ra.push_back(init_ra);
        dec.push_back(init_dec);
        probability.push_back(0.0);
        ranks.push_back(0);
    }

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }

    std::vector<double> probability_sorted;
    std::vector<int> ranks_sorted;
    std::vector<int> indices;
    int start_idx = probability.size() - 1;
    sortNodes(probability, probability_sorted, ranks, ranks_sorted, indices, start_idx);
    chooseTopRegion(costs, probability_sorted, ranks_sorted, indices, costs_filtered, probability_filtered, ranks_filtered, target_sum_prob);

    start_idx = 0;  // after choosing top region
    return start_idx;
}


int buildGraph_force_initpos(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow, int init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }
    return init_pos;
}

int buildGraph_force_initpos(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               int init_pos) {
    std::vector<double> ra, dec;
    read_data_from_file(filename, probability, ranks, ra, dec);

    if (is_deepslow) {
        compute_costs_deepslow(ra, dec, costs, slew_rate, dwell_time);
    }
    else {
        compute_costs_geodesic(ra, dec, costs, slew_rate, dwell_time);
    }

    std::vector<double> probability_sorted;
    std::vector<int> ranks_sorted;
    std::vector<int> indices;
    int start_idx = init_pos;
    sortNodes(probability, probability_sorted, ranks, ranks_sorted, indices, start_idx);
    chooseTopRegion(costs, probability_sorted, ranks_sorted, indices, costs_filtered, probability_filtered, ranks_filtered, target_sum_prob);

    start_idx = 0;  // after choosing top region
    return start_idx;
}



std::vector<int> recoverOriginalIdx(const std::vector<int>& indices, const std::vector<int>& path) {
    // Convert path to original indices
    std::vector<int> original_path;
    for (int node : path) {
        original_path.push_back(indices[node]);
    }
    return original_path;
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

struct Node {
    double x, y;
};

void readOPlibFile(const std::string& filename, 
                   std::vector<std::vector<double>>& costs,
                   std::vector<double>& prizes,
                   std::vector<int>& ranks) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    std::string line;
    int dimension = 0;
    bool reading_coords = false;
    bool reading_scores = false;
    std::vector<Node> nodes;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line.find("DIMENSION") != std::string::npos) {
            auto pos = line.find(":");
            if (pos != std::string::npos) {
                std::string num_str = line.substr(pos + 1);
                dimension = std::stoi(num_str);
                nodes.resize(dimension);
                prizes.resize(dimension, 0.0);
            }
        }
        else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
            reading_coords = true;
            reading_scores = false;
        }
        else if (line.find("NODE_SCORE_SECTION") != std::string::npos) {
            reading_coords = false;
            reading_scores = true;
        }
        else if (line.find("DEPOT_SECTION") != std::string::npos) {
            break; 
        }
        else if (reading_coords) {
            std::istringstream iss(line);
            int idx;
            double x, y;
            iss >> idx >> x >> y;
            if (idx >= 1 && idx <= dimension) {
                nodes[idx - 1] = {x, y};  // 0-index
            }
        }
        else if (reading_scores) {
            std::istringstream iss(line);
            int idx;
            double prize;
            iss >> idx >> prize;
            if (idx >= 1 && idx <= dimension) {
                prizes[idx - 1] = prize;
            }
        }
    }

    infile.close();
    std::cout << "[DEBUG] Reading file: " << filename << std::endl;
    std::cout << "[DEBUG] dimension = " << dimension << std::endl;
    std::cout << "[DEBUG] After reading, prizes.size() = " << prizes.size() << std::endl;


    costs.assign(dimension, std::vector<double>(dimension, 0.0));
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            if (i == j) continue;
            double dx = nodes[i].x - nodes[j].x;
            double dy = nodes[i].y - nodes[j].y;
            costs[i][j] = std::sqrt(dx * dx + dy * dy);
        }
    }

    ranks.resize(dimension);
    std::iota(ranks.begin(), ranks.end(), 0);
}

