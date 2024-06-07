#include "delta_stepping_parallel.h"
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include <cmath>
#include <atomic>
#include <algorithm>

std::vector<double> deltaSteppingParallel(const Graph& graph, int source, double delta, int numThreads) {
    int numVertices = graph.getNumVertices();
    std::vector<std::atomic<double>> distances(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        distances[i] = std::numeric_limits<double>::max();
    }
    distances[source] = 0.0;

    std::vector<std::vector<int>> buckets(numVertices + 1);
    buckets[0].push_back(source);

    std::vector<std::vector<vwPair>> adjacencyList = graph.getAdjacencyList();

    for (int currentBucket = 0; currentBucket <= numVertices; ++currentBucket) {
        if (!buckets[currentBucket].empty()) {
            std::vector<std::thread> threads;
            std::vector<std::vector<int>> localBuckets(numThreads, std::vector<int>());

            // Define thread function within loop to capture local scope variables
            auto processBucket = [&](int startIdx, int endIdx, int threadIndex) {
                for (int idx = startIdx; idx < endIdx; ++idx) {
                    int u = buckets[currentBucket][idx];
                    for (const auto& edge : adjacencyList[u]) {
                        int v = edge.first;
                        double weight = edge.second;
                        double newDist = distances[u] + weight;

                        while (true) {
                            double oldDist = distances[v].load(std::memory_order_relaxed);
                            if (newDist >= oldDist || distances[v].compare_exchange_weak(oldDist, newDist, std::memory_order_relaxed)) {
                                break;
                            }
                        }

                        int bucketIndex = static_cast<int>(newDist / delta);
                        if (bucketIndex > numVertices) bucketIndex = numVertices;
                        localBuckets[threadIndex].push_back(v);
                    }
                }
            };

            // Divide work among threads
            int bucketSize = buckets[currentBucket].size();
            int chunkSize = (bucketSize + numThreads - 1) / numThreads;
            for (int i = 0; i < numThreads; ++i) {
                int startIdx = i * chunkSize;
                int endIdx = std::min(startIdx + chunkSize, bucketSize);
                threads.emplace_back(processBucket, startIdx, endIdx, i);
            }

            // Join threads
            for (auto& thread : threads) {
                thread.join();
            }

            // Merge thread-local buckets into global buckets
            for (int i = 0; i < numThreads; ++i) {
                for (int j = 0; j <= numVertices; ++j) {
                    buckets[j].insert(buckets[j].end(), localBuckets[i].begin(), localBuckets[i].end());
                }
            }
        }
    }

    // Convert atomic<double> to double
    std::vector<double> finalDistances(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        finalDistances[i] = distances[i].load();
    }

    return finalDistances;
}