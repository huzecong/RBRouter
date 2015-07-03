//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_SEQUENTIALEMBEDDING_H
#define RBROUTER_SEQUENTIALEMBEDDING_H


#include <vector>
#include <utility>
#include "RBBase.h"
#include "RBShortestPath.h"


struct RBSENet {
	const int type;
};

struct RBSEAttachedNet : public RBSENet {
	const int type = 0;
};

struct RBSEIncidentNet : public RBSENet {
	const int type = 1;
};

struct RBSEPort {
	double start_angle, end_angle;
};

struct RBSERegion {
	int n_ports, n_nets;
	RBSEPort port[2];
	RBSENet *net[2];
};

class RBSEVertex {
	std::vector<RBSEVertex *> path;
	LinkedList<RBSEAttachedNet *> attached_list;
	CyclicLinkedList<RBSEAttachedNet *> incident_list;
	std::vector<RBSERegion> region_list;
};

class RBSequentialEmbedding {
	typedef std::vector<RBSEVertex>::iterator Point_Iterator;
	RBNet net;
	std::vector<RBSEVertex> point;

	RBShortestPath graph;

public:
	RBSequentialEmbedding() {}
	RBSequentialEmbedding(const RBNet &net,
						  const std::vector<unsigned int> &seq);

	Point_Iterator point_iterator_begin() {
		return point.begin();
	}
	Point_Iterator point_iterator_end() {
		return point.end();
	}

	void roar(const Point_Iterator &it);

	double length();
};


#endif //RBROUTER_SEQUENTIALEMBEDDING_H
