
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <limits>
#include <algorithm>

//simple class for test optimize in future
class Graph {
public:
    int V;
    std::vector<std::vector<double>> costs;

    Graph(int V);
    Graph(const std::string &fileName);
};