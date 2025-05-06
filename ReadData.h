#pragma once

#include <vector>
#include <string>
#include <limits>

//pass init_pos in form {ra, dec}
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow,
               std::pair<double, double> init_pos);

//pass init_pos in form index
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow,
               int init_pos);

//pass init_pos in form {ra, dec}; filtere region with sum prob of target_sum_prob
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               std::pair<double, double> init_pos);

//pass init_pos in form index; filtere region with sum prob of target_sum_prob            
int buildGraph(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               int init_pos);

int buildGraph_force_initpos(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<int>& ranks, std::vector<int>& ranks_filtered,
               double slew_rate, double dwell_time, bool is_deepslow, double target_sum_prob,
               int init_pos);

int buildGraph_force_initpos(const std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<int>& ranks,
               double slew_rate, double dwell_time, bool is_deepslow, int init_pos);

std::vector<int> recoverOriginalIdx(const std::vector<int>& indices, const std::vector<int>& path);

void readOPlibFile(const std::string& filename, 
                   std::vector<std::vector<double>>& costs,
                   std::vector<double>& prizes, std::vector<int>& ranks);

int removeNodes(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping);