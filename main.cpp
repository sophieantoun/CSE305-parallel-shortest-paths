// #include <iostream>
// #include "delta_stepping.h"
// #include <algorithm>
// #include "dijkstra.h"
// #include "graph.h"
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
#include "delta_stepping_parallel.h"


// int main() {
//     Graph graph(6);

//     // Add edges
//     graph.addEdge(0, 1, 50);
//     graph.addEdge(0, 2, 50);
//     graph.addEdge(1, 3, 40);
//     graph.addEdge(2, 4, 70);
//     graph.addEdge(3, 5, 70);
//     graph.addEdge(4, 5, 40);

//     // Create DeltaSteppingSequential object with delta value
//     DeltaSteppingSequential deltaStepping(graph, 0, true);

//     // Compute shortest paths
//     deltaStepping.solveLightHeavy();

//     // Print solution
//     deltaStepping.printSolution();

//     return 0;
// }
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



// typedef std::pair<int, int> Edge; // (vertex, weight)
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
// Function to parse a graph from a file
Graph loadGraphFromFile(int numVertices, const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Failed to open the file." << std::endl;
        exit(1);
    }
    Graph graph(numVertices);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int u, v;
        double weight;
        if (iss >> u >> v >> weight) {
            graph.addEdge(u, v, weight);
        } else {
            std::cerr << "Error: Failed to parse line: " << line << std::endl;
        }
    }
    file.close();
    return graph;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0]
                  << " <N: number of nodes in input graph> "
                  << " <path to input graph file or type [random, grid]> "
                  << " <algorithm to use: [dijkstra, deltastepping, paralleldeltastepping]> "
                  << " <delta for deltastepping algorithm: [float]> "
                  << " <optional: number of edges for random graph> "
                  << " <optional: min weight for random graph> "
                  << " <optional: max weight for random graph> " << std::endl;
        return 1;
    }

    int numVertices = std::stoi(argv[1]);
    std::string graphSource = argv[2];
    std::string algorithm = argv[3];
    double delta = argc > 4 ? std::stod(argv[4]) : 1.0;
    int numEdges = argc > 5 ? std::stoi(argv[5]) : numVertices * 5;
    double minWeight = argc > 6 ? std::stod(argv[6]) : 1.0;
    double maxWeight = argc > 7 ? std::stod(argv[7]) : 10.0;

    Graph graph(numVertices);
    if (graphSource == "random") {
        graph = createRandomGraph(numVertices, numEdges, minWeight, maxWeight);
    } else if (graphSource == "grid") {
        graph = createGridGraph(numVertices, 0.5);
    } else {
        graph = loadGraphFromFile(numVertices, graphSource);
    }

    std::vector<double> distances;

    auto startTime = std::chrono::high_resolution_clock::now();

    if (algorithm == "dijkstra") {
        distances = dijkstra(graph, 0);
    } else if (algorithm == "deltastepping") {
        DeltaSteppingSequential deltaStepping(graph, 0, delta, false);
        deltaStepping.run();
        distances = deltaStepping.getDistances();
    } else if (algorithm == "paralleldeltastepping") {
        distances = deltaSteppingParallel(graph, 0, delta, 2);
    } else {
        std::cerr << "Unsupported algorithm. Use dijkstra, deltastepping, or paralleldeltastepping." << std::endl;
        return 1;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << algorithm << " took " << duration << " ms." << std::endl;

    if (numVertices <= 20) {
        std::cout << "Vertex\tDistance from Source" << std::endl;
        for (int i = 0; i < numVertices; ++i) {
            std::cout << i << "\t" << distances[i] << std::endl;
        }
    }

    std::ofstream outFile("output_" + algorithm + ".txt");
    for (const auto& dist : distances) {
        outFile << dist << "\n";
    }

    // Run Dijkstra for comparison and print results
    std::cout << "Running Dijkstra's algorithm for comparison..." << std::endl;
    startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> dijkstraDistances = dijkstra(graph, 0);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    std::cout << "Dijkstra's algorithm took " << duration << " ms." << std::endl;

    if (numVertices <= 20) {
        std::cout << "Vertex\tDijkstra Distance from Source" << std::endl;
        for (int i = 0; i < numVertices; ++i) {
            std::cout << i << "\t" << dijkstraDistances[i] << std::endl;
        }
    }

    // Compare the results
    bool isCorrect = true;
    for (int i = 0; i < numVertices; ++i) {
        if (distances[i] != dijkstraDistances[i]) {
            std::cout << "Mismatch at vertex " << i << ": Delta Stepping = " << distances[i] << ", Dijkstra = " << dijkstraDistances[i] << std::endl;
            isCorrect = false;
        }
    }

    if (isCorrect) {
        std::cout << "Delta Stepping results match Dijkstra's results!" << std::endl;
    } else {
        std::cout << "Delta Stepping results do not match Dijkstra's results!" << std::endl;
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
