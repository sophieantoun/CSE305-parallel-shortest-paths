#ifndef DELTASTEPPINGPARALLEL2_H
#define DELTASTEPPINGPARALLEL2_H

#include "graph.h"
#include <vector>
#include <set>
#include <mutex>
#include <thread>
#include <iostream>
#include <cmath>
#include "ThreadPool.h"
class DeltaSteppingParallel2 {
public:
    DeltaSteppingParallel2(const Graph& graph, int source, double delta, bool debug, int numThreads);
    void run();

    const std::vector<double> &getDistances() const;

private:
    struct Bucket {
        std::set<int> vertices;
        std::mutex mutex; 
    };
    const Graph& graph;
    int source;
    double delta;
    bool debug;
    int numThreads;
    std::vector<double> distances;
    std::vector<int> predecessors;
    std::vector<std::vector<int>> lightEdges;
    std::vector<std::vector<int>> heavyEdges;
    int bucketIndex;
    std::mutex bucketMutex;
    std::vector<Bucket> buckets;
    ThreadPool pool; 


    void collectEdges(std::vector<std::pair<int, double>>& edges, bool isLight);
    void processBucket(int threadId);
    void classifyEdges();
    void processEdges(bool isLight);
    void relaxSingleEdge(const GraphEdge& edge);  // Correctly handle a single GraphEdge
    void relaxEdges(const std::vector<std::pair<int, double>>& edges);  // Handle multiple edges

    void printBuckets() const;
};

#endif // DELTASTEPPINGPARALLEL2_H
