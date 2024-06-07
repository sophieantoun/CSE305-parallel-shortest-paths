#ifndef DELTA_STEPPING_PARALLEL_H
#define DELTA_STEPPING_PARALLEL_H

#include "graph.h" // Ensure the Graph class is included
#include <list>
#include <limits>
#include <vector>

#include <thread>
#include <mutex>

#include <cmath>
// Declares the deltaSteppingParallel function which calculates shortest paths in a graph using the Delta-Stepping algorithm in a parallel manner.
std::vector<double> deltaSteppingParallel(const Graph& graph, int source, double delta, int numThreads);

#endif