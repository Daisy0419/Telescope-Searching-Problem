//
//  ReadFile.cpp
//  
//
//  Created by Alice Paul on 4/10/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
//
//
#include <cmath>
#include <numeric> 
#include "ReadFile.hpp"

void tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}



bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
    str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool has_prefix(const std::string &str, const std::string &prefix)
{
    auto res = std::mismatch(prefix.begin(), prefix.end(), str.begin());
    return res.first == prefix.end();
}

// bool is_number(const std::string& s)
// {
//     return !s.empty() && std::find_if(s.begin(),
//         s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
// }

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int graphFromFile(const std::string &fileName, Graph &G, std::string &graphName, double &meanEdgeWeight, int &numNodes) {
    std::string line;
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return -1;
    }

    // Read metadata
    std::string ignoredValue, type;
    std::getline(infile, line); // NAME
    std::istringstream iss(line);
    iss >> ignoredValue >> graphName;

    infile >> ignoredValue >> type; // TYPE
    int dimension = 0;
    infile >> ignoredValue >> dimension; // DIMENSION

    // Validate type
    if (type != "AaronTSP") {
        std::cerr << "Error: Unsupported graph type '" << type << "'" << std::endl;
        return -1;
    }

    // Initialize graph vertices based on DIMENSION
    for (int i = 0; i < dimension; ++i) {
        G.addVertex(i);
    }
    numNodes = dimension;

    // Skip to EDGE_WEIGHT_SECTION
    while (std::getline(infile, line)) {
        if (line == "EDGE_WEIGHT_SECTION") {
            break;
        }
    }

    // Parse edges
    double totalEdgeWeight = 0;
    int edgesAdded = 0;
    while (std::getline(infile, line)) {
        std::istringstream edgeStream(line);
        int head, tail;
        double distance;
        if (!(edgeStream >> head >> tail >> distance)) {
            std::cerr << "Error: Malformed edge line '" << line << "'" << std::endl;
            continue;
        }

        G.addEdge(head, tail, distance);
        totalEdgeWeight += distance;
        edgesAdded++;
    }

    if (edgesAdded > 0) {
        meanEdgeWeight = totalEdgeWeight / edgesAdded;
    } else {
        std::cerr << "Warning: No edges added to the graph" << std::endl;
        meanEdgeWeight = 0.0;
    }

    std::cout << "Graph loaded: " << graphName << std::endl
              << "Nodes added: " << numNodes << std::endl
              << "Edges added: " << edgesAdded << std::endl
              << "Mean edge weight: " << meanEdgeWeight << std::endl;

    return 0;
}

// read distances from csv file
std::map<std::pair<int,int>,double> readDistances(const std::string &fileName)
{
    std::map<std::pair<int,int>,double> dists;
    std::ifstream infile;
    infile.open(fileName);
    std::string nextline;
    
    if (infile.is_open())
    {
        while ( getline (infile,nextline) ){
            // Get line and break by spaces
            std::vector<std::string> words;
            tokenize(nextline, words, ",");

            int i = 100*std::stoi(words[0]);
            int j = 100*std::stoi(words[1]);
            double weight = std::stod(words[2]);
            dists[std::make_pair(i,j)] = weight;
            dists[std::make_pair(j,i)] = weight;
        }
        infile.close();
    }
    
    return dists;
}

void computeDistances(const std::vector<std::vector<double>>& costs,
                        std::map<std::pair<int,int>, double>& dists) {
    int dimension = costs.size();
    
    for(int i = 0; i < dimension; ++i) {
        for(int j = i+1; j < dimension; ++j) {
            dists[{i, j}] = costs[i][j];
        }
    }
}


int expandDemandsAsChains( std::map<std::pair<int,int>, double>& dists,
    const std::vector<double>& prizes, int factor) {
    std::vector<int> int_prizes(prizes.size());
    for(int i = 0; i < prizes.size(); ++i) {
        int_prizes[i] = int(prizes[i] * factor);
    }

    int next_id = int_prizes.size();
    for (size_t v = 0; v < prizes.size(); ++v) {
        int v0 = v;
        int prize = int_prizes[v];

        int prev = v0;
        for (int i = 1; i < prize; ++i) {
            int current = next_id++;
            dists[{prev, current}] = 0.0;
            // dists[{current, prev}] = 0.0;
            prev = current;
        }
    }
    return next_id;
}


void expandStartPositon(
    std::map<std::pair<int,int>, double>& dists,
    int start_pos, int& next_id) {

    int v0 = start_pos;
    int prize = next_id;

    int prev = v0;
    for (int i = 1; i < prize; ++i) {
        int current = next_id++;
        dists[{prev, current}] = 0.0;
        // dists[{current, prev}] = 0.0;
        prev = current;
    }
}


int buildGraph(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes, 
                Graph &G, double &meanEdgeWeight, int &numNodes, int init_pos) {
    std::vector<int> prizes_original;
    std::map<std::pair<int,int>, double> distances;
    computeDistances(costs, distances);
    int next_id = expandDemandsAsChains(distances, prizes, 1000);
    expandStartPositon(distances, init_pos, next_id);

    for (int i = 0; i < next_id-1; ++i) {
        G.addVertex(i);
    }

    // Parse edges
    double totalEdgeWeight = 0;
    int edgesAdded = 0;
    for(auto dist : distances) {
        int head=dist.first.first, tail=dist.first.second;
        double distance=dist.second;

        G.addEdge(head, tail, distance);
        totalEdgeWeight += distance;
        edgesAdded++;
    }

    if (edgesAdded > 0) {
        meanEdgeWeight = totalEdgeWeight / edgesAdded;
    } else {
        std::cerr << "Warning: No edges added to the graph" << std::endl;
        meanEdgeWeight = 0.0;
    }

    std::cout << "Nodes added: " << numNodes << std::endl
              << "Edges added: " << edgesAdded << std::endl
              << "Mean edge weight: " << meanEdgeWeight << std::endl;

    return 0;
}
