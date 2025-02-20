#include <gurobi_c++.h>
#include <vector>

struct ProblemData {
    int N;            
    int startIndex;      
    double Budget;
    std::vector<std::vector<double>> Cost; 
    std::vector<double> Prize;        

    ProblemData(int n, int s, double b,
                const std::vector<std::vector<double>>& c,
                const std::vector<double>& p)
      : N(n), startIndex(s), Budget(b), Cost(c), Prize(p)
    {}
};

std::vector<int> solveWithGurobi(const ProblemData& data);
std::vector<int> gurobiSolve(const std::vector<std::vector<double>> &Cost, const std::vector<double> &Prize,
                int start, double Budget);