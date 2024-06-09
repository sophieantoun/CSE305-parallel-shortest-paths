#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <iostream>
#include <set>
#include <utility>
#include <limits>
#include "Edge.h"
#include "/opt/homebrew/Cellar/libomp/18.1.5/include/omp.h"
#include "/Users/sca/opt/anaconda3/include/omp.h"

typedef std::pair<int, double> vwPair; // (vertex, weight)

class Graph
{
private:
    int numVertices;                           // number of vertices
    int numEdges;                          // number of edges
    std::vector<std::vector<vwPair> > adjacencyList; // adjacency list

public:
    // Initialize the graph with a given number of vertices
    Graph(int n) : numVertices(n), numEdges(0)
    {
        adjacencyList.resize(n);
    }

    // Get the number of vertices in the graph
    int getNumVertices() const
    {
        return numVertices;
    }

    // Get the number of edges in the graph
    int getNumEdges() const
    {
        return numEdges;
    }

    // Get the adjacency list of the graph
    std::vector<std::vector<vwPair> > getAdjacencyList() const
    {
        return adjacencyList;
    }

    // Add an edge to the graph between a source and a destination vertex with a given weight
    void addEdge(int source, int destination, double weight)
    {
    
            if (source < 0 || source >= numVertices || destination < 0 || destination >= numVertices)
        {
            std::cerr << "Error: Vertex out of bounds (source: " << source << ", destination: " << destination << ")" << std::endl;
            return;
        }

        adjacencyList[source].push_back(std::make_pair(destination, weight));
        // adjacencyList[destination].push_back(std::make_pair(source, weight));  // uncomment to make the graph undirected
        numEdges++;
    }

    // Get all edges of the graph
    std::vector<Edge> getEdges() const
    {
        std::vector<Edge> edges;
        for (int u = 0; u < numVertices; ++u)
        {
            for (const auto& v : adjacencyList[u])
            {
                edges.emplace_back(u, v.first, v.second);
            }
        }
        return edges;
    }

    // Get the weight of an edge between two vertices
    double getEdgeWeight(int u, int v) const
    {
        for (const auto& edge : adjacencyList[u])
        {
            if (edge.first == v)
            {
                return edge.second;
            }
        }
        return std::numeric_limits<double>::infinity();
    }

    // Get the neighbors of a vertex
    std::set<int> getVertexNeighbors(int u) const
    {
        std::set<int> neighbors;
        for (const auto& edge : adjacencyList[u])
        {
            neighbors.insert(edge.first);
        }
        return neighbors;
    }

    // Compute the maximum degree of the graph
    int getMaxDegree() const
    {
        int maxDegree = 0;
        for (const auto& neighbors : adjacencyList)
        {
            if (neighbors.size() > maxDegree)
            {
                maxDegree = neighbors.size();
            }
        }
        return maxDegree;
    }
};

#endif // GRAPH_H
