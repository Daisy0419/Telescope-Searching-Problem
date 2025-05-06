#pragma once

#include <vector>

//random number generater
double random_double(double min, double max);
int random_int(int min, int max);
std::vector<int> unique_random_ints(int min, int max, int n);

//2-opt
void fix_cross(std::vector<int>& path, const std::vector<std::vector<double>>& costs);
