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
#include "DeltaSteppingParallel2.h"
// #include "/opt/homebrew/Cellar/libomp/18.1.5/include/omp.h"
// #include "/Users/sca/opt/anaconda3/include/omp.h"


// bool compareDistances(const std::vector<double> &distances, const std::vector<double> &expected) {
//     return distances.size() == expected.size() && std::equal(distances.begin(), distances.end(), expected.begin());
// }

// // Function to test correctness
// void testCorrectness() {
//     int passedTests = 0;
//     double delta = 3.0; // Delta for Delta-Stepping

//     // Test case 1: Simple graph with 5 vertices
//     Graph g1(5);
//     g1.addEdge(0, 1, 10);
//     g1.addEdge(0, 3, 5);
//     g1.addEdge(1, 2, 1);
//     g1.addEdge(3, 1, 3);
//     g1.addEdge(3, 4, 2);
//     g1.addEdge(4, 2, 9);
//     g1.addEdge(4, 0, 7);

//     std::vector<double> distances1D = dijkstra(g1, 0);
//     std::vector<double> distances1DS = deltaStepping(g1, 0, delta);

//     std::vector<double> expected1 = {0, 8, 9, 5, 7};
//     if (compareDistances(distances1D, expected1) && compareDistances(distances1DS, expected1)) {
//         ++passedTests;
//     }
    
//     // Test case 2: Graph with negative weights (shouldn't be used with Dijkstra, but to see behavior)
//     Graph g2(4);
//     g2.addEdge(0, 1, 2);
//     g2.addEdge(1, 2, -5); // Negative weight
//     g2.addEdge(2, 3, 1);
//     g2.addEdge(3, 0, 3);

//     std::vector<double> distances2D = dijkstra(g2, 0);
//     std::vector<double> distances2DS = deltaStepping(g2, 0, delta);

//     std::vector<double> expected2 = {0, 2, -3, -2}; // Expected outcome may not be reliable due to negative weights
//     if (compareDistances(distances2D, expected2) && compareDistances(distances2DS, expected2)) {
//         ++passedTests;
//     }
    
//     // Test case 3: Disconnected graph
//     Graph g3(6);
//     g3.addEdge(0, 1, 4);
//     g3.addEdge(1, 2, 6);
//     g3.addEdge(2, 0, 5);
//     // Components 3, 4, 5 are disconnected from 0, 1, 2

//     std::vector<double> distances3D = dijkstra(g3, 0);
//     std::vector<double> distances3DS = deltaStepping(g3, 0, delta);

//     std::vector<double> expected3 = {0, 4, 10, std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
//     if (compareDistances(distances3D, expected3) && compareDistances(distances3DS, expected3)) {
//         ++passedTests;
//     }

//     std::cout << "Number of test cases passed: " << passedTests << "/3" << std::endl;
// }

// int main() {
//     std::cout << "Testing correctness of Dijkstra's and Delta-Stepping algorithms:" << std::endl;
//     testCorrectness();
//     return 0;
// }



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

    std::vector<double> dijkstraDistances, deltaSteppingDistances, parallelDeltaSteppingDistances1, parallelDeltaSteppingDistances2;

    // Print graph details
    std::cout << "Running on graph: " << graphType << " with " << numVertices << " vertices and weights in range [" << minWeight << ", " << maxWeight << "]" << std::endl;

    // Benchmark Dijkstra's algorithm
    auto startTime = std::chrono::high_resolution_clock::now();
    std::cout << "Running Dijkstra's algorithm" << std::endl;
    dijkstraDistances = dijkstra(graph, 0);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Dijkstra's algorithm took " << duration << " ms." << std::endl;

    // Benchmark Delta Stepping Sequential algorithm
    startTime = std::chrono::high_resolution_clock::now();
    DeltaSteppingSequential deltaStepping(graph, 0, delta, false);
    deltaStepping.run();
    deltaSteppingDistances = deltaStepping.getDistances();
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Delta Stepping Sequential algorithm took " << duration << " ms." << std::endl;

    // Create output file
    std::ofstream outFile("output_results.txt", std::ios_base::app);

    // Benchmark Parallel Delta Stepping algorithm with different thread counts
    for (int threads = 2; threads <= maxNumThreads; threads += 2) {
        std::cout << "Running Parallel Delta Stepping algorithm with " << threads << " threads" << std::endl;

        // // Version 1
        startTime = std::chrono::high_resolution_clock::now();
        DeltaSteppingParallel deltaSteppingParallel(graph, 0, delta, false, threads);
        deltaSteppingParallel.run();
        parallelDeltaSteppingDistances1 = deltaSteppingParallel.getDistances();
        endTime = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "Parallel Delta Stepping algorithm (version 1) with " << threads << " threads took " << duration << " ms." << std::endl;

        // Write benchmark details to file
        outFile << "Graph Type: " << graphType << "\n";
        outFile << "Vertices: " << numVertices << "\n";
        outFile << "Edges: " << numEdges << "\n";
        outFile << "Min Weight: " << minWeight << "\n";
        outFile << "Max Weight: " << maxWeight << "\n";
        outFile << "Delta: " << delta << "\n";
        outFile << "Threads: " << threads << "\n";
        outFile << "Algorithm: Parallel Delta Stepping (version 1)\n";
        outFile << "Time: " << duration << " ms\n";

        // Version 2
        startTime = std::chrono::high_resolution_clock::now();
        DeltaSteppingParallel2 deltaSteppingParallel2(graph, 0, delta, false, threads);
        deltaSteppingParallel2.run();
        parallelDeltaSteppingDistances2 = deltaSteppingParallel2.getDistances();
        endTime = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "Parallel Delta Stepping algorithm (version 2) with " << threads << " threads took " << duration << " ms." << std::endl;

        // Write benchmark details to file for version 2
        outFile << "Graph Type: " << graphType << "\n";
        outFile << "Vertices: " << numVertices << "\n";
        outFile << "Edges: " << numEdges << "\n";
        outFile << "Min Weight: " << minWeight << "\n";
        outFile << "Max Weight: " << maxWeight << "\n";
        outFile << "Delta: " << delta << "\n";
        outFile << "Threads: " << threads << "\n";
        outFile << "Algorithm: Parallel Delta Stepping (version 2)\n";
        outFile << "Time: " << duration << " ms\n";

        // Compare the results
        bool isCorrect = true;
        for (int i = 0; i < numVertices; ++i) {
            if (deltaSteppingDistances[i] != dijkstraDistances[i]) {
                outFile << "Mismatch at vertex " << i << ": Delta Stepping = " << deltaSteppingDistances[i] << ", Dijkstra = " << dijkstraDistances[i] << "\n";
                isCorrect = false;
            }
            if (parallelDeltaSteppingDistances1[i] != dijkstraDistances[i]) {
                outFile << "Mismatch at vertex " << i << ": Parallel Delta Stepping 1 = " << parallelDeltaSteppingDistances1[i] << ", Dijkstra = " << dijkstraDistances[i] << "\n";
                isCorrect = false;
            }
            if (parallelDeltaSteppingDistances2[i] != dijkstraDistances[i]) {
                outFile << "Mismatch at vertex " << i << ": Parallel Delta Stepping 2 = " << parallelDeltaSteppingDistances2[i] << ", Dijkstra = " << dijkstraDistances[i] << "\n";
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
    // Define parameters for benchmarks
    std::vector<int> nodeCounts = {200, 500};
    std::vector<std::string> graphTypes = {"dense", "random"};
    std::vector<std::pair<double, double>> weightRanges = {{1, 10}};
    double delta = 3.0; // Example delta value

    for (int numVertices : nodeCounts) {
        for (const std::string& graphType : graphTypes) {
            for (const auto& weightRange : weightRanges) {
                double minWeight = weightRange.first;
                double maxWeight = weightRange.second;
                int numEdges = numVertices * 5;
                runBenchmarks(numVertices, graphType, delta, 8, numEdges, minWeight, maxWeight);
                
            }
        }
    }

    return 0;
}

// // // Graph g(4);
// // // g.addEdge(0, 1, 1);
// // // g.addEdge(1, 2, 1);
// // // g.addEdge(2, 3, 1);
// // // g.addEdge(0, 3, 10); // Longer edge that should not be the shortest path

// // // // Assuming getAdjacencyList returns vector<vector<pair<int, double>>>
// // // auto adjList = g.getAdjacencyList();
// // // for (int i = 0; i < adjList.size(); ++i) {
// // //     std::cout << "Node " << i << " has edges to:\n";
// // //     for (const auto& edge : adjList[i]) {
// // //         std::cout << "  - Node " << edge.first << " with weight " << edge.second << "\n";
// // //     }
// // // } 
// // }
