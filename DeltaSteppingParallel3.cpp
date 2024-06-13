#include "DeltaSteppingParallel3.h"
#include "ThreadPool.h"

DeltaSteppingParallel3::DeltaSteppingParallel3(const Graph& graph, int source, double delta, bool debug, int numThreads)
    : graph(graph), source(source), delta(delta), debug(debug), numThreads(numThreads),
      distances(graph.getNumVertices(), std::numeric_limits<double>::infinity()),
      predecessors(graph.getNumVertices(), -1),
      lightEdges(graph.getNumVertices()), heavyEdges(graph.getNumVertices()),
      bucketIndex(0), buckets(static_cast<int>(graph.getNumVertices() / delta) + 1), 
      pool(numThreads)
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

void DeltaSteppingParallel3::classifyEdges() {
    for (const GraphEdge& edge : graph.getEdges()) {
        if (edge.getWeight() <= delta) {
            lightEdges[edge.getSource()].push_back(edge.getDestination());
        } else {
            heavyEdges[edge.getSource()].push_back(edge.getDestination());
        }
    }
}

void DeltaSteppingParallel3::run() {
    while (bucketIndex < buckets.size()) {
        std::cout << "Current bucket index: " << bucketIndex << ", Size: " << buckets[bucketIndex].vertices.size() << std::endl;
        processEdges(true);  
        processEdges(false); 
        ++bucketIndex;
        while (bucketIndex < buckets.size() && buckets[bucketIndex].vertices.empty()) {
            ++bucketIndex;
        }
    }
}


void DeltaSteppingParallel3::processEdges(bool isLight) {
    std::vector<std::pair<int, double>> edges;
    {
        std::lock_guard<std::mutex> lock(bucketMutex); 
        collectEdges(edges, isLight);
    }

    int threads = std::min(numThreads, static_cast<int>(edges.size()));
    if (threads == 0) threads = 1;  

    int blockSize = edges.size() / threads;
    for (int i = 0; i < threads; ++i) {
        int start = i * blockSize;
        int end = (i + 1 < threads) ? (i + 1) * blockSize : edges.size();
        pool.enqueue([this, start, end, edges]() { // Note: Passing edges by value to avoid access to shared data
            this->relaxEdges(std::vector<std::pair<int, double>>(edges.begin() + start, edges.begin() + end));
        });
    }
}


void DeltaSteppingParallel3::collectEdges(std::vector<std::pair<int, double>>& edges, bool isLight) {
    std::lock_guard<std::mutex> lock(bucketMutex);
    for (int vertex : buckets[bucketIndex].vertices) {
        const auto& edgeList = isLight ? lightEdges[vertex] : heavyEdges[vertex];
        for (int neighbor : edgeList) {
            double weight = graph.getEdgeWeight(vertex, neighbor);
            edges.emplace_back(neighbor, distances[vertex] + weight);
        }
    }
}


void DeltaSteppingParallel3::relaxSingleEdge(const GraphEdge& edge) {
    int fromVertex = edge.getSource();
    int toVertex = edge.getDestination();
    double weight = edge.getWeight();
    double newDistance = distances[fromVertex] + weight;

    if (newDistance < distances[toVertex]) {
        std::lock_guard<std::mutex> lock(bucketMutex);
        int oldBucketIndex = static_cast<int>(std::floor(distances[toVertex] / delta));
        int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

        if (oldBucketIndex >= 0 && oldBucketIndex < buckets.size()) {
            buckets[oldBucketIndex].vertices.erase(toVertex);
        }
        if (newBucketIndex >= 0 && newBucketIndex < buckets.size()) {
            buckets[newBucketIndex].vertices.insert(toVertex);
        }

        distances[toVertex] = newDistance;
        predecessors[toVertex] = fromVertex;

        if (debug) {
            std::cout << "Relaxed edge (" << fromVertex << " -> " << toVertex << ") with new distance " << newDistance << std::endl;
        }
    }
}



void DeltaSteppingParallel3::relaxEdges(const std::vector<std::pair<int, double>>& edges) {
    std::lock_guard<std::mutex> lock(bucketMutex);
    for (const auto& edge : edges) {
        int vertex = edge.first;
        double newDistance = edge.second;

        if (newDistance < distances[vertex]) {
            int oldBucketIndex = static_cast<int>(std::floor(distances[vertex] / delta));
            int newBucketIndex = static_cast<int>(std::floor(newDistance / delta));

            if (oldBucketIndex >= 0 && oldBucketIndex < buckets.size()) {
                buckets[oldBucketIndex].vertices.erase(vertex);
            }
            if (newBucketIndex >= 0 && newBucketIndex < buckets.size()) {
                buckets[newBucketIndex].vertices.insert(vertex);
            }

            distances[vertex] = newDistance;
            predecessors[vertex] = vertex; // Update if needed, otherwise ensure this is correct

            if (debug) {
                std::cout << "Relaxed edge to vertex " << vertex << " with new distance " << newDistance << std::endl;
            }
        }
    }
}



void DeltaSteppingParallel3::processBucket(int threadId) {
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
            relaxSingleEdge(edge);  
        }

        {
            std::lock_guard<std::mutex> lock(buckets[bucketIndex].mutex);
            currentBucket.assign(buckets[bucketIndex].vertices.begin(), buckets[bucketIndex].vertices.end());
        }

        for (const GraphEdge& edge : heavyRequests) {
            relaxSingleEdge(edge);  
        }
    }
}



const std::vector<double>& DeltaSteppingParallel3::getDistances() const {
    return distances;
}

void DeltaSteppingParallel3::printBuckets() const {
    for (size_t bucketId = 0; bucketId < buckets.size(); ++bucketId) {
        std::cout << "Bucket [" << bucketId << "], size " << buckets[bucketId].vertices.size() << ": ";
        for (int vertex : buckets[bucketId].vertices) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}
