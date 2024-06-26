#include "DeltaSteppingParallel.h"

DeltaSteppingParallel::DeltaSteppingParallel(const Graph& graph, int source, double delta, bool debug, int numThreads)
    : graph(graph), source(source), delta(delta), debug(debug), numThreads(numThreads),
      distances(graph.getNumVertices(), std::numeric_limits<double>::infinity()),
      predecessors(graph.getNumVertices(), -1),
      lightEdges(graph.getNumVertices()), heavyEdges(graph.getNumVertices()),
      buckets(static_cast<int>(graph.getNumVertices() / delta) + 1), bucketIndex(0)
{
    classifyEdges();

    for (int i = 0; i < graph.getNumVertices(); ++i) {
        if (i != source)
            buckets.back().insert(i);
    }

    buckets[0].insert(source);
    distances[source] = 0;
    predecessors[source] = source;

}

void DeltaSteppingParallel::classifyEdges() {
    for (const GraphEdge& edge : graph.getEdges()) {
        if (edge.getWeight() <= delta) {
            lightEdges[edge.getSource()].push_back(edge.getDestination());
        } else {
            heavyEdges[edge.getSource()].push_back(edge.getDestination());
        }
    }
}

void DeltaSteppingParallel::run() {
    while (bucketIndex < buckets.size()) {
        std::vector<std::thread> threads;
        for (int i = 0; i < numThreads; ++i) {
            threads.emplace_back(&DeltaSteppingParallel::processBucket, this, i);
        }

        for (auto& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        ++bucketIndex;
        while (bucketIndex < buckets.size() && buckets[bucketIndex].empty()) {
            ++bucketIndex;
        }
    }
}

void DeltaSteppingParallel::processBucket(int threadId) {
    std::set<int> currentBucket;
    {
        std::lock_guard<std::mutex> lock(bucketMutex);
        currentBucket = buckets[bucketIndex];
    }

    while (!currentBucket.empty()) {
        std::vector<GraphEdge> lightRequests;
        std::vector<GraphEdge> heavyRequests;
        for (int vertex : currentBucket) {
            {
                std::lock_guard<std::mutex> lock(bucketMutex);
                buckets[bucketIndex].erase(vertex);
            }

            for (int lightEdge : lightEdges[vertex]) {
                lightRequests.emplace_back(vertex, lightEdge, graph.getEdgeWeight(vertex, lightEdge));
            }

            for (int heavyEdge : heavyEdges[vertex]) {
                heavyRequests.emplace_back(vertex, heavyEdge, graph.getEdgeWeight(vertex, heavyEdge));
            }
        }

        for (const GraphEdge& edge : lightRequests) {
            relaxEdge(edge);
        }

        {
            std::lock_guard<std::mutex> lock(bucketMutex);
            currentBucket = buckets[bucketIndex];
        }
    

        for (const GraphEdge& edge : heavyRequests) {
            relaxEdge(edge);
        }
}
}
void DeltaSteppingParallel::relaxEdge(const GraphEdge& edge) {
    int fromVertex = edge.getSource();
    int toVertex = edge.getDestination();
    double weight = edge.getWeight();
    double newDistance = distances[fromVertex] + weight;

    if (newDistance < distances[toVertex]) {
        int oldBucketIndex = static_cast<int>(std::floor(distances[toVertex] / delta));
        int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

        {
            std::lock_guard<std::mutex> lock(bucketMutex);
            if (oldBucketIndex < buckets.size() && oldBucketIndex >= 0) {
                buckets[oldBucketIndex].erase(toVertex);
            }
            if (newBucketIndex < buckets.size() && newBucketIndex >= 0) {
                buckets[newBucketIndex].insert(toVertex);
            }
        }

        distances[toVertex] = newDistance;
        predecessors[toVertex] = fromVertex;

        if (debug) {
            std::cout << "Relaxed edge (" << fromVertex << " -> " << toVertex << ") with new distance " << newDistance << std::endl;
        }
    }
}

const std::vector<double>& DeltaSteppingParallel::getDistances() const {
    return distances;
}
