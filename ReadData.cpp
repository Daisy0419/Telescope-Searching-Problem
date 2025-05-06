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

    probability.clear();
    ranks.clear();
    ra.clear();
    dec.clear();

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

// edge cost = geodesic distancing + dwell time
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

// edge cost = max axes distancing + dwell time
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

