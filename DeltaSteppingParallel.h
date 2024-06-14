#ifndef DELTASTEPPINGPARALLEL_H
#define DELTASTEPPINGPARALLEL_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include <thread>
#include <atomic>
#include <mutex>
#include "graph.h"

class DeltaSteppingParallel {
public:
    DeltaSteppingParallel(const Graph& graph, int source, double delta, bool debug, int numThreads);

    void run();
    const std::vector<double>& getDistances() const;

private:
    void classifyEdges();
    void processBucket(int threadId);
    void relaxEdge(const GraphEdge& edge);

    const Graph& graph;
    int source;
    double delta;
    bool debug;
    int numThreads;

    std::vector<double> distances;
    std::vector<int> predecessors;
    std::vector<std::vector<int>> lightEdges;
    std::vector<std::vector<int>> heavyEdges;
    std::vector<std::set<int>> buckets;
    int bucketIndex;

    std::mutex bucketMutex;
};

#endif // DELTASTEPPINGPARALLEL_H