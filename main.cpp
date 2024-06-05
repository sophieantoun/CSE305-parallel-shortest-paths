#include <iostream>
#include "delta_stepping.h"
#include <algorithm>
#include "dijkstra.h"
#include "graph.h"

// // Function to compare distances
// bool compareDistances(const std::vector<double> &distances, const std::vector<double> &expected) {
//     return distances.size() == expected.size() && std::equal(distances.begin(), distances.end(), expected.begin());
// }

// // Function to test correctness
// void testCorrectness() {
//     int passedTests = 0;

//     // Test case 1: Simple graph with 5 vertices
//     Graph g1(5);
//     g1.addEdge(0, 1, 10);
//     g1.addEdge(0, 3, 5);
//     g1.addEdge(1, 2, 1);
//     g1.addEdge(3, 1, 3);
//     g1.addEdge(3, 4, 2);
//     g1.addEdge(4, 2, 9);
//     g1.addEdge(4, 0, 7);

//     std::vector<double> distances1 = dijkstra(g1, 0);
//     std::vector<double> expected1 = {0, 8, 9, 5, 7};
//     if (compareDistances(distances1, expected1)) {
//         ++passedTests;
//     }

//     // Test case 2: Graph with negative weights (shouldn't be used with Dijkstra, but to see behavior)
//     Graph g2(4);
//     g2.addEdge(0, 1, 2);
//     g2.addEdge(1, 2, -5);
//     g2.addEdge(2, 3, 1);
//     g2.addEdge(3, 0, 3);

//     std::vector<double> distances2 = dijkstra(g2, 0);
//     std::vector<double> expected2 = {0, 2, -3, -2};
//     if (compareDistances(distances2, expected2)) {
//         ++passedTests;
//     }

//     // Test case 3: Disconnected graph
//     Graph g3(6);
//     g3.addEdge(0, 1, 4);
//     g3.addEdge(1, 2, 6);
//     g3.addEdge(2, 0, 5);
//     g3.addEdge(3, 4, 7);
//     g3.addEdge(4, 5, 3);

//     std::vector<double> distances3 = dijkstra(g3, 0);
//     std::vector<double> expected3 = {0, 4, 10, std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
//     if (compareDistances(distances3, expected3)) {
//         ++passedTests;
//     }

//     std::cout << "Number of test cases passed: " << passedTests << "/3" << std::endl;
// }

// // Function to test Delta-Stepping algorithm
// void testDeltaStepping() {
//     Graph g(5);
//     g.addEdge(0, 1, 10);
//     g.addEdge(0, 3, 5);
//     g.addEdge(1, 2, 1);
//     g.addEdge(3, 1, 3);
//     g.addEdge(3, 4, 2);
//     g.addEdge(4, 2, 9);
//     g.addEdge(4, 0, 7);

//     double delta = 3.0;
//     std::vector<double> distances = deltaStepping(g, 0, delta);

//     for (int i = 0; i < distances.size(); ++i) {
//         if (distances[i] == std::numeric_limits<double>::max()) {
//             std::cout << "Vertex " << i << " is unreachable" << std::endl;
//         } else {
//             std::cout << "Distance to vertex " << i << " is " << distances[i] << std::endl;
//         }
//     }
// }

// int main() {
//     std::cout << "Testing correctness of Dijkstra's algorithm:" << std::endl;
//     testCorrectness();

//     std::cout << "Testing Delta-Stepping algorithm:" << std::endl;
//     testDeltaStepping();

//     return 0;
// }

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

typedef std::pair<int, int> Edge; // (vertex, weight)
const int INF = std::numeric_limits<int>::max();

// Function to generate a random graph
Graph createRandomGraph(int numVertices, int numEdges) {
    Graph graph(numVertices);
    std::srand(std::time(nullptr)); // Seed for random generator
    for (int i = 0; i < numEdges; ++i) {
        int u = std::rand() % numVertices;
        int v = std::rand() % numVertices;
        double weight = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
        graph.addEdge(u, v, weight);
    }
    return graph;
}

// Function to find the smallest and largest edge weights
Edge findMinMaxEdgeWeights(const std::vector<std::vector<vwPair>>& adjList, int numVertices) {
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
    Graph graph(gridSize);
    double weight;
    double probability;
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int node = i * gridSize + j;

            // Right node
            if (j < gridSize - 1) {
                int rightNode = node + 1;
                weight = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                probability = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                if (probability > edgeProbability) {
                    graph.addEdge(node, rightNode, weight);
                }
            }

            // Left node
            if (j > 0) {
                int leftNode = node - 1;
                weight = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                probability = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                if (probability > edgeProbability) {
                    graph.addEdge(node, leftNode, weight);
                }
            }

            // North node
            if (i > 0) {
                int northNode = node - gridSize;
                weight = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                probability = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                if (probability > edgeProbability) {
                    graph.addEdge(node, northNode, weight);
                }
            }

            // South node
            if (i < gridSize - 1) {
                int southNode = node + gridSize;
                weight = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                probability = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // Random double in [0,1]
                if (probability > edgeProbability) {
                    graph.addEdge(node, southNode, weight);
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

    if (argc < 2) {
        std::cout << "Usage: " << argv[0]
                  << " <N: number of nodes in input graph> "
                  << " <path to input graph file:>"
                  << " <algorithm to use: [dijkstra, deltastepping]> "
                  << " <optional: delta for deltastepping algorithm: [float]> " << std::endl;
        return 1;
    }

    // Parse command-line arguments
    int numVertices = std::stoi(argv[1]);
    std::string graphFile = argv[2];
    std::string algorithm = argv[3];

    Graph graph(numVertices);
    if (graphFile == "random") {
        graph = createRandomGraph(numVertices, 150); // We can modify number of edges as needed
    } else if (graphFile == "grid") {
        graph = createGridGraph(numVertices, 0.7); // Modify edge probability as needed
    } else {
        graph = loadGraphFromFile(numVertices, graphFile);
    }

    std::vector<double> dijkstraDistances;

    // Execute the specified algorithm
    std::cout << "Input graph has " << numVertices << " nodes and " << graph.getNumEdges() << " edges" << std::endl;
    if (algorithm == "dijkstra") {
        // Measure execution time for Dijkstra's algorithm
        auto startTime = std::chrono::high_resolution_clock::now();

        dijkstraDistances = dijkstra(graph, 0);
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        std::cout << "Dijkstra's algorithm took " << duration << " µs." << std::endl;

        if(numVertices < 20) {
            std::cout << "Vertex\tDistance from Source" << std::endl;
            for (int i = 0; i < graph.getNumVertices(); ++i) {
                std::cout << i << "\t" << dijkstraDistances[i] << std::endl;
            }
        }

        std::ofstream outFile("output_dijkstra.txt");
        for (const auto &dist : dijkstraDistances) outFile << dist << "\n";
    }

    return 0;
}
