//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_RBMATH_H
#define RBROUTER_RBMATH_H


#include <climits>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>
#include "Box2D/b2BlockAllocator.h"

typedef unsigned int ID;

const double eps = 1e-8;
const double INFI = 1e20;

#define __RB_HAS_ENSURE 1
#define __RB_DEBUG 1

#if __RB_DEBUG
#include <iostream>
#define debug(x) std::cerr << x << std::endl
#else
#define debug(x) ;
#endif

#if __RB_HAS_ENSURE
#define ensure(cond, msg) { if (!(cond)) perror(msg); }
#else
#define ensure(cond, msg) ;
#endif

struct Point {
	double x, y;
	Point() {}
	Point(double _x, double _y) : x(_x), y(_y) {}
};

inline double sqr(double x) {
	return x * x;
}

inline double dist(const Point &a, const Point &b) {
	return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
}

inline double angle(const Point &a, const Point &b) {
	return atan2(b.y - a.y, b.x - a.x);
}

inline bool equal(const double a, const double b) {
	return fabs(a - b) < eps;
}

inline bool lt(const double a, const double b) {
	return b > a + eps;
}

inline bool gt(const double a, const double b) {
	return a > b + eps;
}

inline double cross_product(const Point &o, const Point &a, const Point &b) {
	return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

inline bool colinear(const Point &a, const Point &b, const Point &c) {
	return equal(cross_product(a, b, c), 0.0);
}

inline bool between(const Point &o, const Point &a, const Point &b) {
	return std::min(a.x, b.x) <= o.x && std::max(a.x, b.x) >= o.x
		   && std::min(a.y, b.y) <= o.y && std::max(a.y, b.y) >= o.y;
}

inline int sign(double x) {
	return x > eps ? 1 : x < -eps ? -1 : 0;
}

inline bool intersect_strict(const Point &a, const Point &b,
							 const Point &s, const Point &t) {
	if (sign(cross_product(s, a, t)) * sign(cross_product(s, b, t)) < 0
		&& sign(cross_product(a, s, b)) * sign(cross_product(a, t, b)) < 0) {
		return true;
	} else return false;
}

template<typename T>
inline void print(const std::vector<T> &vec) {
	for (const auto &x : vec)
		std::cerr << x << " ";
	std::cerr << std::endl;
}


template<typename T>
struct LinkedList {

	static b2BlockAllocator *allocator;

	struct ListNode {
		ListNode *prev, *next;
		T data;
		ListNode() {}
		ListNode(const T &_data) : data(_data) {}
	} *head, *tail;
	size_t n_nodes;
	LinkedList() {
		if (this->allocator == NULL) {
			this->allocator = new b2BlockAllocator();
		}

		void *mem = this->allocator->Allocate(sizeof(ListNode));
		this->tail = new (mem) ListNode();
		this->tail->prev = this->tail->next = NULL;
		this->head = this->tail;
		this->n_nodes = 0;
	}
	LinkedList(const LinkedList &list) {
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		this->tail = new (mem) ListNode();
		this->tail->prev = this->tail->next = NULL;
		this->head = this->tail;
		this->n_nodes = 0;
		for (ListNode *p = list.head; p != list.tail; p = p->next)
			this->append(p->data);
	}
	~LinkedList() {
		this->clear();
		this->tail->~ListNode();
		this->allocator->Free(this->tail, sizeof(ListNode));
	}
	const size_t size() const {
		return n_nodes;
	}
	void reset() {
		this->clear();
		this->tail->~ListNode();
		this->allocator->Free(this->tail, sizeof(ListNode));
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		this->tail = new (mem) ListNode();
		this->tail->prev = this->tail->next = NULL;
		this->head = this->tail;
	}
	void clear() {
		for (ListNode *p = this->head; p != this->tail; p = p->next) {
			p->~ListNode();
			this->allocator->Free(p, sizeof(ListNode));
		}
		this->head = this->tail;
		this->tail->prev = this->tail->next = NULL;
		this->n_nodes = 0;
	}
	// Append data before node x
	ListNode *append(ListNode *x, const T &data) {
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		ListNode *q = new (mem) ListNode(data);
		q->next = x;
		if (x->prev != NULL)
			x->prev->next = q;
		q->prev = x->prev, x->prev = q;
		if (this->head == x)
			this->head = q;
		++this->n_nodes;
		return q;
	}
	// Append data at end of queue
	ListNode *append(const T &data) {
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		ListNode *q = new (mem) ListNode(data);
		q->next = this->tail;
		if (this->tail->prev != NULL)
			this->tail->prev->next = q;
		q->prev = this->tail->prev, this->tail->prev = q;
		if (this->head == this->tail)
			this->head = q;
		++this->n_nodes;
		return q;
	}
	// Remove node x
	void remove(ListNode *q) {
		if (q->prev) q->prev->next = q->next;
		if (q->next) q->next->prev = q->prev;
		if (this->head == q) this->head = this->head->next;
		q->~ListNode();
		this->allocator->Free(q, sizeof(ListNode));
		--this->n_nodes;
	}
	// Remove given node, returns true on success
	bool remove(const T &x) {
		for (ListNode *p = this->head; p != this->tail; p = p->next)
			if (p->data == x) {
				this->remove(p);
				return true;
			}
		assert(false);
		return false;
	}
	// Return first node
	const ListNode *front() const {
		return this->head;
	}
	// Locate internal node using user-provided function
	virtual ListNode *find_if(const std::function<bool(const T &)> &func) {
		for (ListNode *p = this->head; p != this->tail; p = p->next)
			if (func(p->data)) return p;
		return NULL;
	}
	// Find internal node which minimizes function
	virtual ListNode *find_minimum(const std::function<double(const T &)>
								   &func) {
		double min_val = DBL_MAX;
		ListNode *ret_node = NULL;
		for (ListNode *p = this->head; p != this->tail; p = p->next) {
			double cur = func(p->data);
			if (cur < min_val) {
				min_val = cur;
				ret_node = p;
			}
		}
		return ret_node;
	}
	// Map function on every node
	virtual void map(const std::function<void(T &)> &func) {
		for (ListNode *p = this->head; p != this->tail; p = p->next)
			func(p->data);
	}
};

template<typename T>
b2BlockAllocator *LinkedList<T>::allocator;


struct RBRoutingPlan {
	std::vector<std::vector<ID>> path;
	double length, width, height;
};

class RBNet {
	friend class RBSequentialEmbedding;
	friend class RBRouter;

	std::vector<Point> point;
	std::vector<std::pair<ID, ID>> net;
	std::vector<std::vector<std::pair<ID, ID>>> link;
	double m_width, m_height;

public:
	RBNet() {}
	RBNet(double _w, double _h, const std::vector<Point> &vec)
			: m_width(_w), m_height(_h), point(vec),
			  link(vec.size(), std::vector<std::pair<ID, ID>>()) {}
	ID add_point(const Point &p) {
		point.push_back(p);
		link.push_back(std::vector<std::pair<ID, ID>>());
		return static_cast<ID>(point.size()) - 1;
	}
	void add_net(ID id1, ID id2) {
		ID net_id = this->n_nets();
		ensure(id1 < this->n_points() && id2 < this->n_points(),
			   "Invalid net: ID incorrect");
		ensure(id1 != id2, "Invalid net: id1 == id2");
		net.emplace_back(id1, id2);
		link[id1].push_back(std::make_pair(id2, net_id));
		link[id2].push_back(std::make_pair(id1, net_id));
	}

	size_t n_points() const { return this->point.size(); }
	size_t n_nets() const { return this->net.size(); }
	const std::vector<std::pair<ID, ID>> &links_from(ID x) const {
		ensure(x < this->n_points(), "ID incorrect");
		return this->link[x];
	}
	std::vector<std::pair<ID, ID>> &links_from(ID x) {
		ensure(x < this->n_points(), "ID incorrect");
		return this->link[x];
	}
	double &width() { return m_width; }
	const double &width() const { return m_width; }
	double &height() { return m_height; }
	const double &height() const { return m_height; }

	// Create component from given IDs
	RBNet subnet(std::vector<ID> vec) const;
	// Combine another net into this net
	void combine(const RBNet &b);
};

#endif //RBROUTER_RBMATH_H
