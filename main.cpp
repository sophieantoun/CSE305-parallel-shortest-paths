#include "graph.h"
#include "dijkstra.h"
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

bool compareDistances(const std::vector<double> &distances, const std::vector<double> &expected) {
    return distances.size() == expected.size() && std::equal(distances.begin(), distances.end(), expected.begin());
}

void testDijkstra() {
    int passedTests = 0;

    // Test case 1: Simple graph with 5 vertices
    Graph g1(5);
    g1.addEdge(0, 1, 10);
    g1.addEdge(0, 3, 5);
    g1.addEdge(1, 2, 1);
    g1.addEdge(3, 1, 3);
    g1.addEdge(3, 4, 2);
    g1.addEdge(4, 2, 9);
    g1.addEdge(4, 0, 7);

    std::vector<double> distances1 = dijkstra(g1, 0);
    std::vector<double> expected1 = {0, 8, 9, 5, 7};
    if (compareDistances(distances1, expected1)) {
        ++passedTests;
    }

    // Test case 2: Graph with negative weights (shouldn't be used with Dijkstra, but to see behavior)
    Graph g2(4);
    g2.addEdge(0, 1, 2);
    g2.addEdge(1, 2, -5);
    g2.addEdge(2, 3, 1);
    g2.addEdge(3, 0, 3);

    std::vector<double> distances2 = dijkstra(g2, 0);
    std::vector<double> expected2 = {0, 2, -3, -2};
    if (compareDistances(distances2, expected2)) {
        ++passedTests;
    }

    // Test case 3: Disconnected graph
    Graph g3(6);
    g3.addEdge(0, 1, 4);
    g3.addEdge(1, 2, 6);
    g3.addEdge(2, 0, 5);
    g3.addEdge(3, 4, 7);
    g3.addEdge(4, 5, 3);

    std::vector<double> distances3 = dijkstra(g3, 0);
    std::vector<double> expected3 = {0, 4, 10, std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    if (compareDistances(distances3, expected3)) {
        ++passedTests;
    }

    std::cout << "Number of test cases passed: " << passedTests << "/3" << std::endl;
}

int main() {
    testDijkstra();
    return 0;
}