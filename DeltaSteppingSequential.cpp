#include "DeltaSteppingSequential.h"

DeltaSteppingSequential::DeltaSteppingSequential(const Graph& graph, const int source, const double delta, const bool debug)
    : graph(graph), source(source), debug(debug), delta(delta),
      distances(graph.getNumVertices(), std::numeric_limits<double>::infinity()), 
      predecessors(graph.getNumVertices(), -1), 
      lightEdges(graph.getNumVertices()), 
      heavyEdges(graph.getNumVertices()), 
      buckets(static_cast<int>(graph.getNumVertices() / delta) + 1), 
      bucketIndex(0) 
{
    classifyEdges();

    if (debug)
        std::cout << "Number of buckets: " << buckets.size() << std::endl;

    for (int i = 0; i < graph.getNumVertices(); ++i) {
        if (i != source)
            buckets.back().insert(i);
    }

    // Initialize source bucket
    buckets[0].insert(source);
    distances[source] = 0;
    predecessors[source] = source;

    if (debug) {
        printEdges();
        printBuckets();
    }
}

void DeltaSteppingSequential::classifyEdges()
{
    for (const Edge& edge : graph.getEdges()) {
        if (edge.getWeight() <= delta) {
            lightEdges[edge.getFrom()].push_back(edge.getTo());
        } else {
            heavyEdges[edge.getFrom()].push_back(edge.getTo());
        }
    }
}

void DeltaSteppingSequential::collectEdgesFromBucket(const std::set<int>& bucket, std::vector<Edge>& lightEdges, std::vector<Edge>& heavyEdges)
{
    for (int vertex : bucket) {
        buckets[bucketIndex].erase(vertex);
        if (debug) {
            std::cout << "Removed " << vertex << " from bucket " << bucketIndex << std::endl;
            std::cout << "Bucket [" << bucketIndex << "], size " << buckets[bucketIndex].size() << ": ";
            for (int v : buckets[bucketIndex]) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }

        for (int lightEdge : this->lightEdges[vertex]) {
            lightEdges.emplace_back(vertex, lightEdge, graph.getEdgeWeight(vertex, lightEdge));
        }

        for (int heavyEdge : this->heavyEdges[vertex]) {
            heavyEdges.emplace_back(vertex, heavyEdge, graph.getEdgeWeight(vertex, heavyEdge));
        }
    }
}

void DeltaSteppingSequential::relaxEdges(std::vector<Edge>& edges)
{
    for (const Edge& edge : edges) {
        int fromVertex = edge.getFrom();
        int toVertex = edge.getTo();
        double weight = edge.getWeight();
        double newDistance = distances[fromVertex] + weight;

        if (newDistance < distances[toVertex]) {
            int oldBucketIndex = static_cast<int>(std::floor(distances[toVertex] / delta));
            int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

            if (oldBucketIndex < buckets.size() && oldBucketIndex >= 0) {
                buckets[oldBucketIndex].erase(toVertex);
            }
            if (newBucketIndex < buckets.size() && newBucketIndex >= 0) {
                buckets[newBucketIndex].insert(toVertex);
            }

            distances[toVertex] = newDistance;
            predecessors[toVertex] = fromVertex;

            if (debug) {
                std::cout << "Relaxed edge (" << fromVertex << " -> " << toVertex << ") with new distance " << newDistance << std::endl;
            }
        }
    }
}

void DeltaSteppingSequential::run()
{
    while (bucketIndex < buckets.size()) {
        std::vector<Edge> lightRequests, heavyRequests;
        std::set<int> currentBucket = buckets[bucketIndex];
        if (debug) {
            std::cout << "Processing bucket " << bucketIndex << ", size " << currentBucket.size() << std::endl;
        }
        while (!currentBucket.empty()) {
            collectEdgesFromBucket(currentBucket, lightRequests, heavyRequests);
            relaxEdges(lightRequests);
            lightRequests.clear();
            currentBucket = buckets[bucketIndex];
        }
        relaxEdges(heavyRequests);
        heavyRequests.clear();
        buckets[bucketIndex].clear(); // Clear the processed bucket
        bucketIndex++;
        while (bucketIndex < buckets.size() && buckets[bucketIndex].empty()) {
            bucketIndex++;
        }
    }
}

void DeltaSteppingSequential::printSolution() const
{
    std::cout << "Distances from source:" << std::endl;
    for (int i = 0; i < graph.getNumVertices(); ++i) {
        std::cout << "Vertex " << i << ": " << distances[i] << std::endl;
    }
}

void DeltaSteppingSequential::printEdges() const
{
    std::cout << "Light edges:" << std::endl;
    for (int i = 0; i < graph.getNumVertices(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (int lightEdge : lightEdges[i]) {
            std::cout << lightEdge << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Heavy edges:" << std::endl;
    for (int i = 0; i < graph.getNumVertices(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (int heavyEdge : heavyEdges[i]) {
            std::cout << heavyEdge << " ";
        }
        std::cout << std::endl;
    }
}

void DeltaSteppingSequential::printBuckets() const
{
    for (size_t bucketId = 0; bucketId < buckets.size(); ++bucketId) {
        std::cout << "Bucket [" << bucketId << "], size " << buckets[bucketId].size() << ": ";
        for (int vertex : buckets[bucketId]) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}

const std::vector<double>& DeltaSteppingSequential::getDistances() const {
    return distances;
}
