#ifndef DELTASTEPPINGSEQUENTIAL_H
#define DELTASTEPPINGSEQUENTIAL_H

#include "graph.h"
#include "Edge.h"
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <iostream>

class DeltaSteppingSequential {
public:
    DeltaSteppingSequential(const Graph& graph, const int source, const bool is_verbose);
    void solve();
    void solveLightHeavy();
    void printSolution() const;

private:
    const Graph& graph_;
    int source_;
    double delta_;
    bool is_verbose_;

    std::vector<double> dist_;
    std::vector<int> pred_;
    std::vector<std::vector<int>> light_edges_;
    std::vector<std::vector<int>> heavy_edges_;
    std::vector<std::set<int>> buckets_;
    int bucket_counter_ = 0;

    void computeLightAndHeavyEdges();
    void findBucketRequests(const std::set<int>& bucket, std::vector<Edge>* light_requests, std::vector<Edge>* heavy_requests);
    void relax(const Edge& selected_edge);
    void resolveRequests(std::vector<Edge>* requests);
    void printLightAndHeavyEdges() const;
    void printAllBuckets() const;
    void printBucket(size_t bucket_id) const;
};

// Add this function declaration
std::vector<double> deltaStepping(const Graph& graph, int source, double delta);


#endif // DELTASTEPPINGSEQUENTIAL_H
