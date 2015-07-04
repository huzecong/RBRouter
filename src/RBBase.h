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

void ensure(const bool cond, const char *msg) {
	if (!cond) perror(msg);
}

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

const double eps = 1e-8;

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

struct Segment {
	Point s, e;
};


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
	}
	~LinkedList() {
		this->clear();
		this->tail->~ListNode();
		this->allocator->Free(this->tail, sizeof(ListNode));
	}
	const size_t size() const {
		return n_nodes;
	}
	void clear() {
		for (ListNode *p = this->head; p != this->tail; ++p) {
			p->~ListNode();
			this->allocator->Free(p, sizeof(ListNode));
		}
		this->head = this->tail;
		this->tail->prev = this->tail->next = NULL;
	}
	// Append data after node x
	ListNode *append(ListNode *x, const T &data) {
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		ListNode *q = new (mem) ListNode(data);
		q->next = x->next;
		if (q->next != NULL)
			q->next->prev = q;
		q->prev = x, x->next = q;
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
		++this->n_nodes;
		return q;
	}
	// Remove node x
	void remove(ListNode *q) {
		if (q->prev) q->prev->next = q->next;
		if (q->next) q->next->prev = q->prev;
		q->~ListNode();
		this->allocator->Free(q, sizeof(ListNode));
		--this->n_nodes;
	}
	// Remove given node, returns true on success
	bool remove(const T &x) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			if (p->data == x) {
				this->remove(p);
				return true;
			}
		return false;
	}
	// Return first node
	const ListNode *front() const {
		return this->head;
	}
	// Locate internal node using user-provided function
	virtual ListNode *find_if(const std::function<bool(const T &)> &func) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			if (func(p->data)) return p;
		return NULL;
	}
	// Find internal node which minimizes function
	virtual ListNode *find_minimum(const std::function<double(const T &)>
								   &func) {
		double min_val = DBL_MAX;
		ListNode *ret_node = NULL;
		for (ListNode *p = this->head; p != this->tail; ++p) {
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
		for (ListNode *p = this->head; p != this->tail; ++p)
			func(p->data);
	}
};

template<typename T>
b2BlockAllocator *LinkedList<T>::allocator;


struct RBRoutingPlan {
	std::vector<std::vector<ID>> path;
	double length;
};

class RBNet {
	friend class RBSequentialEmbedding;

	std::vector<Point> point;
	std::vector<std::pair<ID, ID>> net;
	std::vector<std::vector<ID>> link;
	double m_width, m_height;

public:
	RBNet() {}
	RBNet(double _w, double _h, const std::vector<Point> &vec)
			: m_width(_w), m_height(_h), point(vec),
			  link(vec.size(), std::vector<ID>()) {}
	ID add_point(const Point &p) {
		point.push_back(p);
		link.push_back(std::vector<ID>());
		return static_cast<ID>(point.size()) - 1;
	}
	void add_net(ID id1, ID id2) {
		ensure(id1 < this->n_points() && id2 < this->n_points(),
			   "Invalid net: ID incorrect");
		ensure(id1 != id2, "Invalid net: id1 == id2");
		net.emplace_back(id1, id2);
		link[id1].push_back(id2);
		link[id2].push_back(id1);
	}

	size_t n_points() const { return this->point.size(); }
	size_t n_nets() const { return this->net.size(); }
	const std::vector<ID> &links_from(ID x) const {
		ensure(x < this->n_points(), "ID incorrect");
		return this->link[x];
	}
	std::vector<ID> &links_from(ID x) {
		ensure(x < this->n_points(), "ID incorrect");
		return this->link[x];
	}
	double &width() { return m_width; }
	const double &width() const { return m_width; }
	double &height() { return m_height; }
	const double &height() const { return m_height; }

	// Create component from given IDs
	RBNet subnet(std::vector<ID> vec) const;
};

#endif //RBROUTER_RBMATH_H
