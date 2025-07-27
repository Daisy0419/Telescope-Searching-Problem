#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>

ILOSTLBEGIN

struct CplexProblemData {
    int N;
    int startIndex;
    int endIndex;
    double Budget;
    std::vector<std::vector<double>> Cost;
    std::vector<double> Prize;

    CplexProblemData(int n, int s, int e, double b,
                const std::vector<std::vector<double>>& c,
                const std::vector<double>& p)
        : N(n), startIndex(s), endIndex(e), Budget(b), Cost(c), Prize(p) {}
};

std::vector<int> solveWithCplexST(const CplexProblemData& data, 
                                double mipGap, double timeLimit,
                                const std::vector<int>& initialRoute,                                
                                double& obj_bound, double& result_gap);


std::vector<int> cplexSolveST(const std::vector<std::vector<double>>& Cost,
                               const std::vector<double>& Prize,
                               int start, int end, double Budget,
                               double& obj_value, double& result_gap,
                               double mipGap, double timeLimit,
                               const std::vector<int>& initialRoute);