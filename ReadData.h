#pragma once

#include "helplers.h"

#include <vector>
#include <string>

double angularDistance(double ra1, double dec1, double ra2, double dec2);
void read_data (std::string& filename, std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, double telescope_speed, double dwell_time);
void read_data_deep_slow(std::string& filename, std::vector<std::vector<double>>& costs, 
               std::vector<double>& probability, double telescope_speed, double dwell_time);

std::vector<std::vector<double>> initialize_distances(int NUM_CITIES);
std::vector<std::vector<double>> initialize_distances(const std::string &fileName);