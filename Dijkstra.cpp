#include <queue>
#include <limits>
#include <stdexcept>

#include "graph.h"

// Dijkstra algorithm 
std::vector<double> dijkstra(const Graph &graph, int source)
{
    int numVertices = graph.getNumVertices();

    // Check that the source vertex is valid
    if (source < 0 || source >= numVertices)
        throw std::invalid_argument("The source vertex is out of the range of the graph");


    std::vector<double> distances(numVertices, std::numeric_limits<double>::max());  //We initialize distances to all vertices as infinite

    std::vector<std::vector<vwPair>> adjacencyList = graph.getAdjacencyList();    //Adjacency list from the graph

    std::priorityQueue<vwPair> priorityQueue; //Priority queue to hold the vertices to be explored

    // Start with the source vertex
    distances[source] = 0.0;
    priorityQueue.push(std::make_pair(0.0, source));

    // Start the Dijkstra algorithm
    while (!priorityQueue.empty())
    {
        // Get the vertex with the smallest distance
        int u = priorityQueue.top().second;
        priorityQueue.pop();

        // Iterate through all the neighboring vertices of u
        for (const auto &edge : adjacencyList[u])
        {
            int v = edge.first;
            double weight = edge.second;

            // Calculate the possible alternative distance to vertex v through u
            double altDistance = distances[u] + weight;

            // Update the distance if it is shorter
            if (altDistance < distances[v])
            {
                distances[v] = altDistance;
                priorityQueue.push(std::make_pair(altDistance, v));
            }
        }
    }

    return distances;
}

