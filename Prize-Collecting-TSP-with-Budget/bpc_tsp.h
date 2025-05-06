#pragma once

#include <vector>

std::vector<int> BPC_TSP(const std::vector<std::vector<double>> &costs, 
            const std::vector<double>& prize, double D, int init_pos);

std::vector<int> BPC_TSP2(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                double D,
                int init_pos);