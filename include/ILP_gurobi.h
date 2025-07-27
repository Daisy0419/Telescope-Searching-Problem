#pragma once

#include <gurobi_c++.h>
#include <vector>
#include <fstream>
#include <iomanip>  
#include <cmath> 

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




class LogToFileCallback : public GRBCallback {
    std::ostream& out_;
    static constexpr long STEP = 5000; 
public:
    explicit LogToFileCallback(std::ostream& o) : out_(o) {
        out_ << "#node,best,bound,gap\n";
    }

protected:
    void callback() override {
        if (where == GRB_CB_MIP) {
            long   node  = static_cast<long>(getDoubleInfo(GRB_CB_MIP_NODCNT));
            if (node % STEP == 0) {  
                double best  = getDoubleInfo(GRB_CB_MIP_OBJBST);
                double bound = getDoubleInfo(GRB_CB_MIP_OBJBND);
                double gap   = 100.0 * (best - bound) /
                               (std::fabs(best) + 1e-9);

                out_ << node << ',' << best << ',' << bound << ','
                     << gap << '\n';
                out_.flush();
            }
        }
        else if (where == GRB_CB_MIPSOL) {
            double obj   = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
            double bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);

            out_ << "#new_incumbent," << obj << ',' << bound << '\n';
            out_.flush();
        }
    }
};
