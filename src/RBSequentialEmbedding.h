//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_SEQUENTIALEMBEDDING_H
#define RBROUTER_SEQUENTIALEMBEDDING_H


#include <vector>
#include <utility>
#include "RBBase.h"
#include "RBShortestPath.h"
#include "Box2D/b2BlockAllocator.h"

enum RBSENetType {
	NullType = -1,
	AttachedNet = 0,
	IncidentNet = 1
};

struct RBSENet {
	const RBSENetType type = NullType;
	ID net_no, from_node;
	RBSENet *counterpart;
	double angle;
};

struct RBSEAttachedNet : public RBSENet {
	const RBSENetType type = AttachedNet;
};

struct RBSEIncidentNet : public RBSENet {
	const RBSENetType type = IncidentNet;
};

struct RBSEPort {
	double start_angle, end_angle;
	RBSEPort() {}
	RBSEPort(double s, double e) : start_angle(s), end_angle(e) {}
	inline bool operator ==(const RBSEPort &p) const {
		return equal(this->start_angle, p.start_angle)
			   && equal(this->end_angle, p.end_angle);
	}
};

struct RBSERegion {
	static b2BlockAllocator *allocator;

	ID id;
	LinkedList<RBSEPort> port;
	LinkedList<RBSENet *> net;
	LinkedList<std::pair<RBSERegion *, ID>> link;
	struct RBSEVertex *belong_vertex;
	bool is_open;

	static ID id_cnt;
	static std::vector<RBSERegion *> regions;
	static std::vector<ID> free_id;

	static void init() {
		id_cnt = 0;
		regions.clear();
		free_id.clear();
	}

	void *operator new(size_t) throw(std::bad_alloc) {
		return allocator->Allocate(sizeof(RBSERegion));
	}

	void operator delete(void *m) {
		allocator->Free(m, sizeof(RBSERegion));
	}

	RBSERegion() {
		if (free_id.size() > 0) {
			id = free_id[free_id.size() - 1];
			free_id.pop_back();
		} else {
			id = id_cnt++;
		}
		if (regions.size() <= id)
			regions.resize(id + 1);
		regions[id] = this;
	}
	~RBSERegion() {
		free_id.push_back(id);
		regions[id] = NULL;
	}
};

struct RBSEVertex {
	static b2BlockAllocator *allocator;

	Point point;
	LinkedList<RBSEAttachedNet *> attached_list;
	LinkedList<RBSEAttachedNet *> incident_list;
	LinkedList<RBSERegion *> region_list;
	LinkedList<RBSERegion *> open_region;

	void *operator new(size_t) throw(std::bad_alloc) {
		return allocator->Allocate(sizeof(RBSEVertex));
	}

	void operator delete(void *m) {
		allocator->Free(m, sizeof(RBSEVertex));
	}
};

class RBSequentialEmbedding {
	typedef std::vector<RBSEVertex *>::iterator Point_Iterator;
	RBNet net;
	std::vector<RBSEVertex *> vertex;

	RBShortestPath graph;

	RBRoutingPlan plan;

public:
	RBSequentialEmbedding() {}
	RBSequentialEmbedding(const RBNet &net,
						  const std::vector<unsigned int> &seq);

	Point_Iterator point_iterator_begin() {
		return vertex.begin();
	}
	Point_Iterator point_iterator_end() {
		return vertex.end();
	}

	void roar(const Point_Iterator &it);

	double length();
};


#endif //RBROUTER_SEQUENTIALEMBEDDING_H
