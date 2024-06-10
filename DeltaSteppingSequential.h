#ifndef DELTASTEPPINGSEQUENTIAL_H
#define DELTASTEPPINGSEQUENTIAL_H

#include <vector>
#include <set>
#include <iostream>
#include <limits>
#include <cmath>
#include "graph.h"
#include "Edge.h"

class DeltaSteppingSequential {
public:
    DeltaSteppingSequential(const Graph& graph, const int source, const double delta, const bool debug);
    void run();
    const std::vector<double>& getDistances() const;
    void printSolution() const;

private:
    const Graph& graph;
    const int source;
    const double delta;
    const bool debug;
    std::vector<double> distances;
    std::vector<int> predecessors;
    std::vector<std::vector<int>> lightEdges;
    std::vector<std::vector<int>> heavyEdges;
    std::vector<std::set<int>> buckets;
    int bucketIndex;

    void initialize();
    void processBucket();
    void processEdgeTasks(std::vector<GraphEdge>& edgeTasks);
    void printEdges() const;
    void printBuckets() const;
};

#endif // DELTASTEPPINGSEQUENTIAL_H
