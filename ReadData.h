#pragma once

#include "helplers.h"

#include <vector>
#include <string>

double angularDistance(double ra1, double dec1, double ra2, double dec2);
void read_data (std::string& filename, std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, double telescope_speed, double dwell_time);

void read_data_filtered(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<double>& ranks, 
               double telescope_speed, double dwell_time, double select_probability);

void read_data_deep_slow(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, double telescope_speed, double dwell_time);

void read_data_deep_slow_filtered(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, std::vector<std::vector<double>>& costs_filtered, 
               std::vector<double>& probability_filtered, std::vector<double>& ranks, 
               double telescope_speed, double dwell_time, double select_probability);

std::vector<std::vector<double>> initialize_distances(int NUM_CITIES);
std::vector<std::vector<double>> initialize_distances(const std::string &fileName);