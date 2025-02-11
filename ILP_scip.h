#pragma once
#include "objscip/objscip.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <memory>

// Problem data 
struct problemData {
    int N;
    int startIndex;
    double Budget;
    std::vector<std::vector<double>> Cost;
    std::vector<double> Prize;

    problemData(int n, int s, double b, const std::vector<std::vector<double>>& cost, const std::vector<double>& prize)
        : N(n), startIndex(s), Budget(b), Cost(cost), Prize(prize) {}
};

class tspProblem : public scip::ObjProbData {
public:
    std::shared_ptr<problemData> data;
    std::vector<std::vector<SCIP_VAR*>> x;
    std::vector<SCIP_VAR*> is_end;
    std::vector<SCIP_VAR*> u;

    tspProblem(std::shared_ptr<problemData> d) : data(d) {}
    virtual SCIP_RETCODE scip_create(SCIP* scip);
};

int solveProblem(const std::vector<std::vector<double>>& Cost,
                 const std::vector<double>& Prize,
                 double Budget, int startIndex);


