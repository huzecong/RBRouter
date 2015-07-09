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
	RBSENet *counterpart;
	double angle;
	RBSENet() {}
	RBSENet(RBSENetType _t) : type(_t) {}
};

struct RBSEAttachedNet : public RBSENet {
	RBSEAttachedNet() : RBSENet(AttachedNet) {}
};

struct RBSEIncidentNet : public RBSENet {
	RBSEIncidentNet() : RBSENet(IncidentNet) {}
};

struct RBSEPort {
	static RBSEPort PortNotFound;
	double start_angle, end_angle;
	RBSEPort() {}
	RBSEPort(double s, double e) : start_angle(s), end_angle(e) {}
	inline bool operator ==(const RBSEPort &p) const {
		return equal(this->start_angle, p.start_angle)
			   && equal(this->end_angle, p.end_angle);
	}
	bool contains_eq(double angle) const;
	bool contains(double angle) const;
	bool contains(const RBSEPort &other) const;
	double length() const;
};

struct RBSERegion {
	static b2BlockAllocator *allocator;

	ID id;
	LinkedList<RBSEPort> port;
	LinkedList<RBSENet *> net;
	LinkedList<std::pair<RBSERegion *, ID>> link;
	struct RBSEVertex *vertex;
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
		: port(region.port), net(region.net), link(region.link),
		  vertex(region.vertex), is_open(region.is_open) {
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
	~RBSERegion() {
		free_id.push_back(this->id);
		regions[this->id] = NULL;
	}
};

struct RBSEVertex {
	static b2BlockAllocator *allocator;

	ID id;
	Point point;
//	LinkedList<RBSEAttachedNet *> attached_list;
//	LinkedList<RBSEAttachedNet *> incident_list;
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

	std::vector<ID> insert_open_edge(RBSEVertex *vertex);

	void roar(ID x);

public:
	RBSequentialEmbedding() {}
	RBSequentialEmbedding(const RBNet &net,
						  const std::vector<unsigned int> &seq);

	const double length() const {
		return plan.length;
	}

	const RBRoutingPlan &routing_plan() const {
		return plan;
	}
};


#endif //RBROUTER_SEQUENTIALEMBEDDING_H
