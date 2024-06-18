#ifndef DELTASTEPPINGPARALLEL2_H
#define DELTASTEPPINGPARALLEL2_H

#include "graph.h"
#include <vector>
#include <set>
#include <mutex>
#include <thread>
#include <iostream>
#include <cmath>


struct Bucket {
    std::set<int> vertices;
    std::mutex mutex;
};

class DeltaSteppingParallel2 {
public:
    DeltaSteppingParallel2(const Graph& graph, int source, double delta, bool debug, int numThreads);
    void run();
    const std::vector<double>& getDistances() const;

private:
    const Graph& graph;
    int source;
    double delta;
    bool debug;
    int numThreads;
    std::vector<double> distances;
    std::vector<int> predecessors;
    std::vector<std::vector<int>> lightEdges;
    std::vector<std::vector<int>> heavyEdges;
    std::vector<Bucket> buckets;
    int bucketIndex;
    std::mutex bucketMutex;

    void classifyEdges();
    void processEdges(bool isLight, std::vector<std::thread>& workers);
    void collectEdges(std::vector<std::pair<int, double>>& edges, bool isLight);
    void relaxEdges(const std::vector<std::pair<int, double>>& edges);
    void processBucket(int threadId);
    void relaxEdge(const GraphEdge& edge);
};

#endif // DELTASTEPPINGPARALLEL2_H