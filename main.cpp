#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "graph.h"
#include "dijkstra.h"
#include "DeltaSteppingSequential.h"
#include "DeltaSteppingParallel.h"
#include "DeltaSteppingParallel3.h"
#include "DeltaSteppingParallel2.h"
#include "ThreadPool.h"

// #include "/opt/homebrew/Cellar/libomp/18.1.5/include/omp.h"
// #include "/Users/sca/opt/anaconda3/include/omp.h"

typedef std::pair<int, int> MinMaxEdge; // (vertex, weight)
const int INF = std::numeric_limits<int>::max();

Graph createRandomGraph(int numVertices, int numEdges, double minWeight, double maxWeight) {
    Graph graph(numVertices);
    std::srand(std::time(nullptr)); // Seed for random generator
    for (int i = 0; i < numEdges; ++i) {
        int u = std::rand() % numVertices;
        int v = std::rand() % numVertices;
        double weight = minWeight + (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (maxWeight - minWeight); // Random double in [minWeight, maxWeight]
        graph.addEdge(u, v, weight);
    }
    return graph;
}

// Function to find the smallest and largest edge weights
MinMaxEdge findMinMaxEdgeWeights(const std::vector<std::vector<vwPair>>& adjList, int numVertices) {
    std::cout << "Finding min and max edge weights" << std::endl;
    int minWeight = INF;
    int maxWeight = -INF;
    for (int i = 0; i < numVertices; ++i) {
        for (const auto& edge : adjList[i]) {
            double weight = edge.second;
            if (weight > maxWeight) {
                maxWeight = weight;
            }
            if (weight < minWeight && weight > 0) {
                minWeight = weight;
            }
        }
    }
    return std::make_pair(minWeight, maxWeight);
}


// Function to generate a grid graph
Graph createGridGraph(int gridSize, double edgeProbability) {
    Graph graph(gridSize * gridSize);
    double weight;
    double probability;
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int node = i * gridSize + j;

            // Right node
            if (j < gridSize - 1) {
                int rightNode = node + 1;
                weight = 1 + (static_cast<double>(rand()) / RAND_MAX) * 9; // Random double in [1, 10]
                probability = static_cast<double>(rand()) / RAND_MAX; // Random double in [0,1]
                if (probability < edgeProbability) {
                    graph.addEdge(node, rightNode, weight);
                    graph.addEdge(rightNode, node, weight); // Ensure bidirectional
                }
            }

            // Left node
            if (j > 0) {
                int leftNode = node - 1;
                weight = 1 + (static_cast<double>(rand()) / RAND_MAX) * 9; // Random double in [1, 10]
                probability = static_cast<double>(rand()) / RAND_MAX; // Random double in [0,1]
                if (probability < edgeProbability) {
                    graph.addEdge(node, leftNode, weight);
                    graph.addEdge(leftNode, node, weight); // Ensure bidirectional
                }
            }

            // North node
            if (i > 0) {
                int northNode = node - gridSize;
                weight = 1 + (static_cast<double>(rand()) / RAND_MAX) * 9; // Random double in [1, 10]
                probability = static_cast<double>(rand()) / RAND_MAX; // Random double in [0,1]
                if (probability < edgeProbability) {
                    graph.addEdge(node, northNode, weight);
                    graph.addEdge(northNode, node, weight); // Ensure bidirectional
                }
            }

            // South node
            if (i < gridSize - 1) {
                int southNode = node + gridSize;
                weight = 1 + (static_cast<double>(rand()) / RAND_MAX) * 9; // Random double in [1, 10]
                probability = static_cast<double>(rand()) / RAND_MAX; // Random double in [0,1]
                if (probability < edgeProbability) {
                    graph.addEdge(node, southNode, weight);
                    graph.addEdge(southNode, node, weight); // Ensure bidirectional
                }
            }
        }
    }
    return graph;
}
Graph createDenseGraph(int numVertices, double minWeight, double maxWeight, double edgeProbability = 0.9) {
    Graph graph(numVertices);
    std::srand(std::time(nullptr)); 
    for (int u = 0; u < numVertices; ++u) {
        for (int v = 0; v < numVertices; ++v) {
            if (u != v && (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) < edgeProbability) {
                double weight = minWeight + (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (maxWeight - minWeight);
                graph.addEdge(u, v, weight);
            }
        }
    }
    return graph;
}
// Function to parse a graph from a file
Graph loadGraphFromFile(int numVertices, const std::string& filename) {
    std::cout << "Loading graph from file: " << filename << std::endl;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Failed to open the file: " << filename << std::endl;
        exit(1);
    }
    Graph graph(numVertices);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int u, v;
        double weight;
        std::cout << "Reading edge: " << line << std::endl; // Debug line
        if (iss >> u >> v >> weight) {
            graph.addEdge(u, v, weight);
        } else {
            std::cerr << "Error: Failed to parse line: " << line << std::endl;
        }
    }
    file.close();
    return graph;
}

        // Version 1
        // startTime = std::chrono::high_resolution_clock::now();
        // DeltaSteppingParallel deltaSteppingParallel(graph, 0, delta, false, threads);
        // deltaSteppingParallel.run();
        // parallelDeltaSteppingDistances1 = deltaSteppingParallel.getDistances();
        // endTime = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // std::cout << "Parallel Delta Stepping algorithm (version 1) with " << threads << " threads took " << duration << " ms." << std::endl;
        // outFile << "Algorithm: Parallel Delta Stepping (version 1)\n";
        // outFile << "Threads: " << threads << "\n";
        // outFile << "Time: " << duration << " ms\n\n";


        // Version 2
        // startTime = std::chrono::high_resolution_clock::now();
        // DeltaSteppingParallel2 deltaSteppingParallel2(graph, 0, delta, false, threads);
        // deltaSteppingParallel2.run();
        // parallelDeltaSteppingDistances2 = deltaStproeppingParallel2.getDistances();
        // endTime = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // std::cout << "Parallel Delta Stepping algorithm (version 2) with " << threads << " threads took " << duration << " ms." << std::endl;
        // outFile << "Algorithm: Parallel Delta Stepping (version 2)\n";
        // outFile << "Threads: " << threads << "\n";
        // outFile << "Time: " << duration << " ms\n\n";


void runBenchmarks(int numVertices, const std::string& graphType, double delta, int maxNumThreads, int numEdges, double minWeight, double maxWeight) {
    Graph graph(numVertices);
    if (graphType == "random") {
        graph = createRandomGraph(numVertices, numEdges, minWeight, maxWeight);
    } else if (graphType == "grid") {
        graph = createGridGraph(numVertices, 0.5);
    } else if (graphType == "dense") {
        graph = createDenseGraph(numVertices, minWeight, maxWeight);
    } else {
        std::cerr << "Unsupported graph type." << std::endl;
        return;
    }

    std::vector<double> dijkstraDistances, deltaSteppingDistances, parallelDeltaSteppingDistances3;

    // Print graph details
    std::cout << "Running on graph: " << graphType << " with " << numVertices << " vertices and weights in range [" << minWeight << ", " << maxWeight << "]" << std::endl;

    // Create output file
    std::ofstream outFile("output_results.txt", std::ios_base::app);

    // Benchmark Dijkstra's algorithm
    auto startTime = std::chrono::high_resolution_clock::now();
    std::cout << "Running Dijkstra's algorithm" << std::endl;
    dijkstraDistances = dijkstra(graph, 0);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Dijkstra's algorithm took " << duration << " ms." << std::endl;
    outFile << "Graph Type: " << graphType << "\n";
    outFile << "Vertices: " << numVertices << "\n";
    outFile << "Edges: " << numEdges << "\n";
    outFile << "Min Weight: " << minWeight << "\n";
    outFile << "Max Weight: " << maxWeight << "\n";
    outFile << "Delta: " << delta << "\n";
    outFile << "Algorithm: Dijkstra\n";
    outFile << "Time: " << duration << " ms\n\n";

    // Benchmark Delta Stepping Sequential algorithm
    startTime = std::chrono::high_resolution_clock::now();
    DeltaSteppingSequential deltaStepping(graph, 0, delta, false);
    deltaStepping.run();
    deltaSteppingDistances = deltaStepping.getDistances();
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Delta Stepping Sequential algorithm took " << duration << " ms." << std::endl;
    outFile << "Algorithm: Delta Stepping Sequential\n";
    outFile << "Time: " << duration << " ms\n\n";

    // Benchmark Parallel Delta Stepping algorithm with different thread counts
    for (int threads = 2; threads <= maxNumThreads; threads += 2) {
        std::cout << "Running Parallel Delta Stepping algorithm with " << threads << " threads" << std::endl;

        // Version 3
        startTime = std::chrono::high_resolution_clock::now();
        std::unique_ptr<DeltaSteppingParallel3> deltaSteppingParallel3(new DeltaSteppingParallel3(graph, 0, delta, false, threads));
        deltaSteppingParallel3->run();
        parallelDeltaSteppingDistances3 = deltaSteppingParallel3->getDistances();
        endTime = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "Parallel Delta Stepping algorithm (version 3) with " << threads << " threads took " << duration << " ms." << std::endl;
        outFile << "Algorithm: Parallel Delta Stepping (version 3)\n";
        outFile << "Threads: " << threads << "\n";
        outFile << "Time: " << duration << " ms\n\n";

        // Compare the results
        bool isCorrect = true;
        for (int i = 0; i < numVertices; ++i) {
            if (deltaSteppingDistances[i] != dijkstraDistances[i]) {
                outFile << "Mismatch at vertex " << i << ": Delta Stepping = " << deltaSteppingDistances[i] << ", Dijkstra = " << dijkstraDistances[i] << "\n";
                isCorrect = false;
            }
            if (parallelDeltaSteppingDistances3[i] != dijkstraDistances[i]) {
                outFile << "Mismatch at vertex " << i << ": Parallel Delta Stepping 3 = " << parallelDeltaSteppingDistances3[i] << ", Dijkstra = " << dijkstraDistances[i] << "\n";
                isCorrect = false;
            }
        }

        if (isCorrect) {
            std::cout << "All algorithms produced matching results!" << std::endl;
            outFile << "All algorithms produced matching results!\n";
        } else {
            outFile << "There were mismatches between the algorithm results!\n";
        }
        outFile << "\n";
    }

    outFile.close();
}

int main() {
    try {
        // Define parameters for benchmarks
        std::vector<int> nodeCounts = {1000, 5000, 10000};
        std::vector<std::string> graphTypes = {"dense", "random", "grid"};
        std::vector<std::pair<double, double>> weightRanges = {{1, 5}, {1, 10}, {1, 100}};
        std::vector<double> deltas = {0.5, 1, 3, 5};
        
        for (int numVertices : nodeCounts) {
            for (const std::string& graphType : graphTypes) {
                for (const auto& weightRange : weightRanges) {
                    for (double delta : deltas) {
                        double minWeight = weightRange.first;
                        double maxWeight = weightRange.second;
                        int numEdges = numVertices * 5;
                        runBenchmarks(numVertices, graphType, delta, 16, numEdges, minWeight, maxWeight);
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "An exception occurred: " << e.what() << std::endl;
    }

    return 0;
}


// void testSimpleGraph() {
//     Graph graph(4);
//     graph.addEdge(0, 1, 1);
//     graph.addEdge(1, 2, 1);
//     graph.addEdge(2, 3, 1);
//     graph.addEdge(0, 3, 10);

//     DeltaSteppingParallel3 deltaStepping(graph, 0, 1.0, true, 4);
//     deltaStepping.run();
//     const std::vector<double>& distances = deltaStepping.getDistances();

//     for (size_t i = 0; i < distances.size(); ++i) {
//         std::cout << "Distance to vertex " << i << " is " << distances[i] << std::endl;
//     }

//     assert(distances[3] == 3); // Check if the distance is as expected
// }


// int main() {
//     testSimpleGraph();
//     return 0;
// }