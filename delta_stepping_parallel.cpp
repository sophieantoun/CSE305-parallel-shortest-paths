#include "delta_stepping_parallel.h"

std::vector<double> deltaSteppingParallel(const Graph& graph, int source, double delta) {
    int numVertices = graph.getNumVertices();
    std::vector<double> distances(numVertices, std::numeric_limits<double>::max());
    std::vector<std::vector<int>> buckets(numVertices + 1);  // Use vector of vectors for better random access
    distances[source] = 0.0;
    buckets[0].push_back(source);

    std::vector<std::vector<vwPair>> adjacencyList = graph.getAdjacencyList();

    int currentBucket = 0;

    while (currentBucket <= numVertices) {
        if (!buckets[currentBucket].empty()) {
            // Parallel processing of each vertex in the current bucket
            #pragma omp parallel
            {
                std::vector<int> localLightEdges;
                std::vector<int> localHeavyEdges;

                #pragma omp for nowait // Distribute vertices to threads
                for (int idx = 0; idx < buckets[currentBucket].size(); ++idx) {
                    int u = buckets[currentBucket][idx];
                    for (const auto& edge : adjacencyList[u]) {
                        int v = edge.first;
                        double weight = edge.second;

                        double newDist = distances[u] + weight;

                        // Critical section to update shared resources
                        #pragma omp critical
                        {
                            if (newDist < distances[v]) {
                                distances[v] = newDist;
                                if (weight <= delta) {
                                    localLightEdges.push_back(v);
                                } else {
                                    localHeavyEdges.push_back(v);
                                }
                            }
                        }
                    }
                }

                // Merge local lists into global lists safely
                #pragma omp critical
                {
                    for (int v : localLightEdges) buckets[static_cast<int>(distances[v] / delta)].push_back(v);
                    for (int v : localHeavyEdges) buckets[static_cast<int>(distances[v] / delta)].push_back(v);
                }
            }
        }
        ++currentBucket;
    }

    return distances;
}