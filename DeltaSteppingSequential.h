#ifndef DELTASTEPPINGSEQUENTIAL_H
#define DELTASTEPPINGSEQUENTIAL_H

#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <iostream>
#include "graph.h"

class DeltaSteppingSequential {
public:
    DeltaSteppingSequential(const Graph& graph, const int source, const double delta, const bool debug = false);

    void run();
    void printSolution() const;
    const std::vector<double>& getDistances() const;

private:
    void classifyEdges();
    void collectEdgesFromBucket(const std::set<int>& bucket, std::vector<Edge>& lightEdges, std::vector<Edge>& heavyEdges);
    void relaxEdges(std::vector<Edge>& edges);
    void printEdges() const;
    void printBuckets() const;

    const Graph& graph;
    int source;
    bool debug;
    double delta;
    std::vector<double> distances;
    std::vector<int> predecessors;
    std::vector<std::vector<int>> lightEdges;
    std::vector<std::vector<int>> heavyEdges;
    std::vector<std::set<int>> buckets;
    int bucketIndex;
};

#endif // DELTASTEPPINGSEQUENTIAL_H
