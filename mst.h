#pragma once
#include <vector>
#include <string>
std::vector<int> mstNaive(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);             

std::vector<int> mstNaiveUpdate(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);                  

std::vector<int> mstLemon3(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);

std::vector<int> mstLemon2(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);

std::vector<int> mstLemon(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);


std::vector<int> mstLemonUpdate(const std::vector<std::vector<double>>& costs, 
                          std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);