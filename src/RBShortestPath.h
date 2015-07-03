//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_RBSHORTESTPATH_H
#define RBROUTER_RBSHORTESTPATH_H


#include <vector>
#include <algorithm>
#include "RBBase.h"

class RBShortestPath {
	typedef unsigned int ID;
	typedef std::pair<ID, double> Edge;

	unsigned int n_vertex, n_edges;

	// Prefer a linked list implementation

public:
	RBShortestPath();
	// Add an edge a->b with cost w, and return the internal ID of the edge
	ID add_edge(ID a, ID b, double w);
	// Remove edge with internal ID x
	void remove_edge(ID x);

	// Find the shortest path from vertex a to b, and return the vertexes on
	// the shortest path
	std::vector<ID> shortest_path(ID a, ID b);
};


#endif //RBROUTER_RBSHORTESTPATH_H
