#pragma once

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

std::vector<int> solveWithGurobi(const ProblemData& data, double mipGap = -1, double timeLimit = -1,
                                 const std::vector<int>& initialRoute = {});


std::vector<int> gurobiSolve(const std::vector<std::vector<double>>& Cost, 
                             const std::vector<double>& Prize,
                             int start, double Budget, 
                             double mipGap = -1, double timeLimit = -1,
                             const std::vector<int>& initialRoute = {});

struct ProblemDataST {
    int N;
    int startIndex;
    int endIndex;
    double Budget;
    std::vector<std::vector<double>> Cost;
    std::vector<double> Prize;

    ProblemDataST(int N_, int start_, int end_, double Budget_,
                const std::vector<std::vector<double>>& Cost_,
                const std::vector<double>& Prize_)
        : N(N_), startIndex(start_), endIndex(end_), Budget(Budget_),
          Cost(Cost_), Prize(Prize_) {}
};

std::vector<int> gurobiSolveST(const std::vector<std::vector<double>>& Cost,
                               const std::vector<double>& Prize,
                               int start, int end, double Budget,
                               double mipGap = -1, double timeLimit = -1,
                               const std::vector<int>& initialRoute = {});