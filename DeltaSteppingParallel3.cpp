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
       // std::cout << "Running bucket index: " << bucketIndex << ", Vertex count: " << buckets[bucketIndex].vertices.size() << std::endl;
        processEdges(true);  
      //  std::cout << "First Process Edge" << std::endl;
        processEdges(false);
        ++bucketIndex;
       // std::cout << "Second" << std::endl;
        while (bucketIndex < buckets.size() && buckets[bucketIndex].vertices.empty()) {
            ++bucketIndex;
       //     std::cout << "Skipping empty bucket index: " << bucketIndex << std::endl;
        }
    }
   // std::cout << "Exiting run loop." << std::endl;
}


void DeltaSteppingParallel3::processEdges(bool isLight) {
    //std::cout << "Processing edges, isLight: " << isLight << ", Bucket Index: " << bucketIndex << std::endl;
    std::vector<std::pair<int, double>> edges;
    {
        // std::lock_guard<std::mutex> lock(bucketMutex);
      //  std::cout << "Collecting edges... Total vertices in bucket: " << buckets[bucketIndex].vertices.size() << std::endl;
        collectEdges(edges, isLight);
    }

    if (edges.empty()) {
        //std::cout << "No edges to process, exiting." << std::endl;
        return;
    }

    int threads = std::min(numThreads, static_cast<int>(edges.size()));
    //std::cout << "Using " << threads << " threads." << std::endl;
    int blockSize = edges.size() / threads;
    std::vector<std::thread> threadsVector;

    for (int i = 0; i < threads; ++i) {
        int start = i * blockSize;
        int end = (i + 1 < threads) ? (i + 1) * blockSize : edges.size();
        threadsVector.emplace_back([this, start, end, &edges]() {
            this->relaxEdges(std::vector<std::pair<int, double>>(edges.begin() + start, edges.begin() + end));
        });
    }

    for (auto& th : threadsVector) {
        th.join();
    }
    //std::cout << "Finished processing edges." << std::endl;
}


void DeltaSteppingParallel3::collectEdges(std::vector<std::pair<int, double>>& edges, bool isLight) {
   //std::cout << "Starting to collect edges. Is light: " << isLight << std::endl;
    std::vector<int> localVertices;

    // Manual locking and unlocking for more control over mutex handling
    //std::cout << "Attempting to lock bucketMutex..." << std::endl;
    bucketMutex.lock();
    //std::cout << "bucketMutex locked." << std::endl;
    localVertices.assign(buckets[bucketIndex].vertices.begin(), buckets[bucketIndex].vertices.end());
    bucketMutex.unlock();
    //std::cout << "bucketMutex unlocked. Collected " << localVertices.size() << " vertices from bucket index " << bucketIndex << std::endl;

    // Immediately check for entrance into vertex loop
    //std::cout << "Checking if we are entering the vertex processing loop..." << std::endl;

    for (int vertex : localVertices) {
       // std::cout << "Processing vertex " << vertex << std::endl;
        std::vector<std::pair<int, double>> tempEdges;
        std::vector<int> edgeList;

        // Lock just before accessing shared edge lists
        bucketMutex.lock();
        edgeList = isLight ? lightEdges[vertex] : heavyEdges[vertex];
        bucketMutex.unlock();

       // std::cout << "Vertex " << vertex << " has " << edgeList.size() << " edges." << std::endl;

        for (int neighbor : edgeList) {
            double weight = graph.getEdgeWeight(vertex, neighbor);
            tempEdges.emplace_back(neighbor, distances[vertex] + weight);
          //  std::cout << "Processed edge from " << vertex << " to " << neighbor << " with weight " << weight << std::endl;
        }

        edges.insert(edges.end(), tempEdges.begin(), tempEdges.end());
    }

//    std::cout << "Finished collecting edges. Total edges collected: " << edges.size() << std::endl;
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


// Example to reduce lock scope and ensure lock release
void DeltaSteppingParallel3::processBucket(int threadId) {
    while (true) {
        std::vector<int> currentBucket;
        {
            std::lock_guard<std::mutex> lock(bucketMutex);
            if (buckets[bucketIndex].vertices.empty()) break;
            currentBucket.assign(buckets[bucketIndex].vertices.begin(), buckets[bucketIndex].vertices.end());
            buckets[bucketIndex].vertices.clear(); // Clear once copied
        }

        std::vector<GraphEdge> lightRequests, heavyRequests;
        for (int vertex : currentBucket) {
            // Assuming lightEdges and heavyEdges are not modified by other threads at this point
            for (int lightEdge : lightEdges[vertex]) {
                lightRequests.emplace_back(vertex, lightEdge, graph.getEdgeWeight(vertex, lightEdge));
            }
            for (int heavyEdge : heavyEdges[vertex]) {
                heavyRequests.emplace_back(vertex, heavyEdge, graph.getEdgeWeight(vertex, heavyEdge));
            }
        }

        // Process requests outside the lock
        for (const GraphEdge& edge : lightRequests) {
            relaxSingleEdge(edge);
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

