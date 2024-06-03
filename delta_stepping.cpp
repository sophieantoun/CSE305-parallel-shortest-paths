#include "delta_stepping.h"

std::vector<double> deltaStepping(const Graph& graph, int source, double delta) {
    int numVertices = graph.getNumVertices();
    std::vector<double> distances(numVertices, std::numeric_limits<double>::max());
    std::vector<std::list<int>> buckets(numVertices + 1); // Buckets for vertices
    distances[source] = 0.0;
    buckets[0].push_back(source);

    std::vector<std::vector<vwPair>> adjacencyList = graph.getAdjacencyList(); // Get a copy of the adjacency list

    int currentBucket = 0;

    while (currentBucket <= numVertices) {
        if (!buckets[currentBucket].empty()) {
            std::list<int> lightEdges;
            std::list<int> heavyEdges;

            // Process each vertex in the current bucket
            while (!buckets[currentBucket].empty()) {
                int u = buckets[currentBucket].front();
                buckets[currentBucket].pop_front();

                for (const auto& edge : adjacencyList[u]) {
                    int v = edge.first;
                    double weight = edge.second;

                    double newDist = distances[u] + weight;

                    if (weight <= delta) {
                        if (newDist < distances[v]) {
                            distances[v] = newDist;
                            lightEdges.push_back(v);
                        }
                    } else {
                        if (newDist < distances[v]) {
                            distances[v] = newDist;
                            heavyEdges.push_back(v);
                        }
                    }
                }
            }

            // Relax light edges
            for (int u : lightEdges) {
                for (const auto& edge : adjacencyList[u]) {
                    int v = edge.first;
                    double weight = edge.second;

                    double newDist = distances[u] + weight;

                    if (newDist < distances[v]) {
                        distances[v] = newDist;
                        int bucketIndex = static_cast<int>(newDist / delta);
                        if (bucketIndex >= numVertices) {
                            bucketIndex = numVertices;
                        }
                        buckets[bucketIndex].push_back(v);
                    }
                }
            }

            // Relax heavy edges
            for (int u : heavyEdges) {
                for (const auto& edge : adjacencyList[u]) {
                    int v = edge.first;
                    double weight = edge.second;

                    double newDist = distances[u] + weight;

                    if (newDist < distances[v]) {
                        distances[v] = newDist;
                        int bucketIndex = static_cast<int>(newDist / delta);
                        if (bucketIndex >= numVertices) {
                            bucketIndex = numVertices;
                        }
                        buckets[bucketIndex].push_back(v);
                    }
                }
            }
        }
        ++currentBucket;
    }

    return distances;
}
