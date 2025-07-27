#pragma once
#include <string> 
#include <vector>

//print path to stdout
double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> rank, std::vector<int> best_path, double dwell_time);

void save_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, double slew_rate, double dwell_time, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> rank,
                  std::vector<int> best_path, 
                  double elapsed_time, double padding);

void save_result2(const std::string& filename, 
                  const std::string& method, 
                  const std::vector<double>& budgets,
                  const std::vector<double>& dwell_time, 
                  const std::vector<std::vector<double>>& costs, 
                  const std::vector<double>& probability, 
                  std::vector<int> rank,
                  std::vector<int> best_path, 
                  double elapsed_time);
                  
// void test_algorithms (std::string file, std::string out_file, double budget, 
//                                 double slew_rate, bool is_deepslow=false);


void test_algorithms_small_instances (std::string file, std::string out_file, 
                    double budget, double slew_rate, bool is_deepslow=false);
void test_algorithms_large_instances (std::string file, std::string out_file, 
                    double budget, double slew_rate, bool is_deepslow=false);
void run_multi_deadlines(std::string file, std::string out_file, std::string out_file2, 
                                std::vector<double> budgets, 
                                double slew_rate, bool is_deepslow=false);