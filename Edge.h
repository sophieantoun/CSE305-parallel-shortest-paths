#ifndef EDGE_H
#define EDGE_H

class GraphEdge {
public:
    GraphEdge(int source, int destination, double weight)
        : sourceVertex(source), destinationVertex(destination), edgeWeight(weight) {}

    int getSource() const { return sourceVertex; }
    int getDestination() const { return destinationVertex; }
    double getWeight() const { return edgeWeight; }

private:
    int sourceVertex;
    int destinationVertex;
    double edgeWeight;
};

#endif // EDGE_H
