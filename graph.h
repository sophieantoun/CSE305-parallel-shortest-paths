#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <iostream>

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
        adjacencyList[source].push_back(std::make_pair(destination, weight));
        // adjacencyList[destination].push_back(std::make_pair(source, weight));  // uncomment to make the graph undirected
        numEdges++;
    }
};

#endif // GRAPH_H
