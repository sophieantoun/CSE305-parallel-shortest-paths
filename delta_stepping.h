#ifndef DELTA_STEPPING_H
#define DELTA_STEPPING_H

#include <vector>
#include <list>
#include <limits>
#include "graph.h"
// Define a pair representing a vertex and weight
typedef std::pair<int, double> vwPair;


// Declaration of the deltaStepping function
std::vector<double> deltaStepping(const Graph& graph, int source, double delta);

#endif // DELTA_STEPPING_H
