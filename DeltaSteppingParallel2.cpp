#include "DeltaSteppingParallel2.h"

DeltaSteppingParallel2::DeltaSteppingParallel2(const Graph& graph, int source, double delta, bool debug, int numThreads)
    : graph(graph), source(source), delta(delta), debug(debug), numThreads(numThreads),
      distances(graph.getNumVertices(), std::numeric_limits<double>::infinity()),
      predecessors(graph.getNumVertices(), -1),
      lightEdges(graph.getNumVertices()), heavyEdges(graph.getNumVertices()),
      buckets(static_cast<int>(graph.getNumVertices() / delta) + 1), bucketIndex(0)
{
    classifyEdges();

    for (int i = 0; i < graph.getNumVertices(); ++i) {
        if (i != source)
            buckets.back().vertices.insert(i);
    }
    
    distances[source] = 0;
    predecessors[source] = source;
    buckets[0].vertices.insert(source);
    

    if (debug) {
        std::cout << "Number of buckets: " << buckets.size() << std::endl;
        printBuckets();
    }
}

void DeltaSteppingParallel2::classifyEdges() {
    for (const GraphEdge& edge : graph.getEdges()) {
        if (edge.getWeight() <= delta) {
            lightEdges[edge.getSource()].push_back(edge.getDestination());
        } else {
            heavyEdges[edge.getSource()].push_back(edge.getDestination());
        }
    }
}

void DeltaSteppingParallel2::run() {
    while (bucketIndex < buckets.size()) {
        std::vector<std::thread> workers;

        processEdges(true, workers);  // Process light edges
        processEdges(false, workers); // Process heavy edges

        // Join all threads
        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }

        // Move to the next non-empty bucket
        ++bucketIndex;
        while (bucketIndex < buckets.size() && buckets[bucketIndex].vertices.empty()) {
            ++bucketIndex;
        }
    }
}

void DeltaSteppingParallel2::processEdges(bool isLight, std::vector<std::thread>& workers) {
    std::vector<std::pair<int, double>> edges;
    collectEdges(edges, isLight);

    int threads = std::min(numThreads, static_cast<int>(edges.size()));
    if (threads == 0) threads = 1;

    int blockSize = edges.size() / threads;
    for (int i = 0; i < threads - 1; ++i) {
        int start = i * blockSize;
        int end = (i + 1) * blockSize;
        workers.emplace_back(&DeltaSteppingParallel2::relaxEdges, this, std::vector<std::pair<int, double>>(edges.begin() + start, edges.begin() + end));
    }
    // Last block
    workers.emplace_back(&DeltaSteppingParallel2::relaxEdges, this, std::vector<std::pair<int, double>>(edges.begin() + (threads - 1) * blockSize, edges.end()));
}

void DeltaSteppingParallel2::collectEdges(std::vector<std::pair<int, double>>& edges, bool isLight) {
    std::lock_guard<std::mutex> lock(bucketMutex);
    for (int vertex : buckets[bucketIndex].vertices) {
        const auto& edgeList = isLight ? lightEdges[vertex] : heavyEdges[vertex];
        for (int neighbor : edgeList) {
            double weight = graph.getEdgeWeight(vertex, neighbor);
            edges.emplace_back(neighbor, distances[vertex] + weight);
        }
    }
}

void DeltaSteppingParallel2::relaxEdges(const std::vector<std::pair<int, double>>& edges) {
    for (const auto& edge : edges) {
        int vertex = edge.first;
        double newDistance = edge.second;

        std::lock_guard<std::mutex> lock(bucketMutex);
        if (newDistance < distances[vertex]) {
            int oldBucketIndex = static_cast<int>(std::floor(distances[vertex] / delta));
            int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

            if (oldBucketIndex < buckets.size() && oldBucketIndex >= 0) {
                buckets[oldBucketIndex].vertices.erase(vertex);
            }
            if (newBucketIndex < buckets.size() && newBucketIndex >= 0) {
                buckets[newBucketIndex].vertices.insert(vertex);
            }

            distances[vertex] = newDistance;
            predecessors[vertex] = source;

            if (debug) {
                std::cout << "Relaxed edge (" << source << " -> " << vertex << ") with new distance " << newDistance << std::endl;
            }
        }
    }
}

void DeltaSteppingParallel2::processBucket(int threadId) {
    std::vector<int> currentBucket;
    {
        std::lock_guard<std::mutex> lock(buckets[bucketIndex].mutex);
        currentBucket.assign(buckets[bucketIndex].vertices.begin(), buckets[bucketIndex].vertices.end());
    }

    while (!currentBucket.empty()) {
        std::vector<GraphEdge> lightRequests;
        std::vector<GraphEdge> heavyRequests;
        
        for (size_t i = 0; i < currentBucket.size(); ++i) {
            int vertex = currentBucket[i];
            {
                std::lock_guard<std::mutex> lock(buckets[bucketIndex].mutex);
                buckets[bucketIndex].vertices.erase(vertex);
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
            std::lock_guard<std::mutex> lock(buckets[bucketIndex].mutex);
            currentBucket.assign(buckets[bucketIndex].vertices.begin(), buckets[bucketIndex].vertices.end());
        }

        for (const GraphEdge& edge : heavyRequests) {
            relaxEdge(edge);
        }
    }
}


void DeltaSteppingParallel2::relaxEdge(const GraphEdge& edge) {
    int fromVertex = edge.getSource();
    int toVertex = edge.getDestination();
    double weight = edge.getWeight();
    double newDistance = distances[fromVertex] + weight;

    if (newDistance < distances[toVertex]) {
        int oldBucketIndex = static_cast<int>(std::floor(distances[toVertex] / delta));
        int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

        {
            std::lock_guard<std::mutex> lock(buckets[bucketIndex].mutex);
            if (oldBucketIndex < buckets.size() && oldBucketIndex >= 0) {
                buckets[oldBucketIndex].vertices.erase(toVertex);
            }
            if (newBucketIndex < buckets.size() && newBucketIndex >= 0) {
                buckets[newBucketIndex].vertices.insert(toVertex);
            }
        }

        distances[toVertex] = newDistance;
        predecessors[toVertex] = fromVertex;

        if (debug) {
            std::cout << "Relaxed edge (" << fromVertex << " -> " << toVertex << ") with new distance " << newDistance << std::endl;
        }
    }
}

const std::vector<double>& DeltaSteppingParallel2::getDistances() const {
    return distances;
}

void DeltaSteppingParallel2::printBuckets() const {
    for (size_t bucketId = 0; bucketId < buckets.size(); ++bucketId) {
        std::cout << "Bucket [" << bucketId << "], size " << buckets[bucketId].vertices.size() << ": ";
        for (int vertex : buckets[bucketId].vertices) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}
