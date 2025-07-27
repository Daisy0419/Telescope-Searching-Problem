#pragma once

#include <vector>
#include <string>
#include <limits>

std::tuple<int, int, double> buildGraphOrienteering(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks, std::vector<double>& dwell_times,
               double slew_rate, bool is_deepslow, int init_pos);

std::pair<int, int> removeNodesOrienteering(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& dwell_times,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping,
                int start_idx, int end_idx);

// void readOPlibFile(const std::string& filename, 
//                    std::vector<std::vector<double>>& costs,
//                    std::vector<double>& prizes, std::vector<int>& ranks);