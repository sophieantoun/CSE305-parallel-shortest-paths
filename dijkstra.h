#ifndef DIJKSTRA_H
#define DIJKSTRA_H
#include <queue>
#include <limits>
#include <stdexcept>
#include <vector>
#include "graph.h"

std::vector<double> dijkstra(const Graph &graph, int source);

#endif // DIJKSTRA_H
