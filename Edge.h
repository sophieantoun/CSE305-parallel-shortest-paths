#ifndef EDGE_H
#define EDGE_H

class Edge {
public:
    Edge(int from, int to, double weight)
        : from_(from), to_(to), weight_(weight) {}

    int getFrom() const { return from_; }
    int getTo() const { return to_; }
    double getWeight() const { return weight_; }

private:
    int from_;
    int to_;
    double weight_;
};

#endif // EDGE_H
