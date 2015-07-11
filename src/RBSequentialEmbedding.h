//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_SEQUENTIALEMBEDDING_H
#define RBROUTER_SEQUENTIALEMBEDDING_H


#include <utility>
#include <vector>
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
	ID net_no;
	RBSENet *counterpart, *opposite;
	struct RBSEVertex *from_vertex, *to_vertex;
	RBSENet() {}
	RBSENet(RBSENetType _t) : type(_t) {
		this->counterpart = this->opposite = NULL;
	}
};

struct RBSEAttachedNet : public RBSENet {
	RBSEAttachedNet() : RBSENet(AttachedNet) {}
};

struct RBSEIncidentNet : public RBSENet {
	RBSEIncidentNet() : RBSENet(IncidentNet) {}
};

struct RBSEPort {
	Point s, e;
	RBSEPort() {}
	RBSEPort(const Point &_s, const Point &_e) : s(_s), e(_e) {}
	bool contains_eq(const Point &p) const;
	bool contains(const Point &p) const;
};

struct RBSERegion {
	static b2BlockAllocator *allocator;

	ID id;
	RBSEAttachedNet *outer_net[2], *inner_net[2];
	RBSEIncidentNet *incident_net[2];
	LinkedList<std::pair<RBSERegion *, ID>> link;
	struct RBSEVertex *vertex;

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
//		if (free_id.size() > 0) {
		if (false) {
			this->id = free_id[free_id.size() - 1];
			free_id.pop_back();
		} else {
			this->id = id_cnt++;
			if (regions.size() <= this->id)
				regions.resize(this->id + 1);
		}
		regions[this->id] = this;
	}
	RBSERegion(const RBSERegion &region)
		: link(region.link) {
//		if (free_id.size() > 0) {
		for (int i = 0; i < 2; ++i) {
			this->inner_net[i] = region.inner_net[i];
			this->outer_net[i] = region.outer_net[i];
			this->incident_net[i] = region.incident_net[i];
		}
		this->vertex = region.vertex;

		if (false) {
			this->id = free_id[free_id.size() - 1];
			free_id.pop_back();
		} else {
			this->id = id_cnt++;
			if (regions.size() <= this->id)
				regions.resize(this->id + 1);
		}
		regions[this->id] = this;
	}
	~RBSERegion() {
		free_id.push_back(this->id);
		regions[this->id] = NULL;
	}
};

struct RBSEVertex {
	static b2BlockAllocator *allocator;

	ID id;
	Point point;
	LinkedList<RBSENet *> attached_list, opposite_list;
	LinkedList<RBSENet *> net_list;
	LinkedList<RBSERegion *> region_list;

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

	std::vector<ID> insert_open_edge(RBSEVertex *vertex);

	void preprocess();
	void roar(ID x);

public:
	RBSequentialEmbedding() {}
	RBSequentialEmbedding(const RBNet &net,
						  const std::vector<unsigned int> &seq);
	~RBSequentialEmbedding();

	const double length() const {
		return plan.length;
	}

	const RBRoutingPlan &routing_plan() const {
		return plan;
	}
};


#endif //RBROUTER_SEQUENTIALEMBEDDING_H
