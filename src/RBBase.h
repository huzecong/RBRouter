//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_RBMATH_H
#define RBROUTER_RBMATH_H


#include <cstdio>
#include <vector>
#include "Box2D/b2BlockAllocator.h"

void ensure(const bool cond, const char *msg) {
	if (!cond) perror(msg);
}

struct Point {
	double x, y;
	Point() {}
	Point(double _x, double _y) : x(_x), y(_y) {}
};

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
	int n_nodes;
	LinkedList() {
		if (this->allocator == NULL) {
			this->allocator = new b2BlockAllocator();
		}

		void *mem = this->allocator->Allocate(sizeof(ListNode));
		this->tail = new (mem) ListNode();
		this->tail->prev = this->tail->next = NULL;
		this->head = this->tail;
	}
	// Append data after node x
	ListNode *append(ListNode *x, const T &data) {
		void *mem = this->allocator->Allocate(sizeof(ListNode));
		ListNode *q = new (mem) ListNode(data);
		q->next = x->next;
		if (q->next != NULL)
			q->next->prev = q;
		q->prev = x, x->next = q;
		++n_nodes;
		return q;
	}
	// Remove node x
	void remove(ListNode *q) {
		if (q->prev) q->prev->next = q->next;
		if (q->next) q->next->prev = q->prev;
		q->~ListNode();
		this->allocator->Free(q, sizeof(ListNode));
		--n_nodes;
	}
	// Locate internal node using user-provided function
	ListNode *find_if(bool (*func)(const T &)) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			if (func(p->data)) return p;
		return NULL;
	}
	// Map function on every node
	void map(void (*func)(T &)) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			func(p->data);
	}
};

template<typename T>
struct CyclicLinkedList : public LinkedList<T> {
	typedef typename LinkedList<T>::ListNode ListNode;

	CyclicLinkedList() : LinkedList<T>() {
		this->head->prev = this->head->next = this->head;
	}
	// Locate internal node using user-provided function
	ListNode *find_if(bool (*func)(const T &)) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			if (func(p->data)) return p;
		if (func(this->tail)) return this->tail;
		return NULL;
	}
	// Map function on every node
	void map(void (*func)(T &)) {
		for (ListNode *p = this->head; p != this->tail; ++p)
			func(p->data);
		func(this->tail);
	}
};


class RoutingPlan {

};

class RBNet {
public:
	typedef unsigned int ID;

private:
	std::vector<Point> point;
	std::vector<std::pair<ID, ID>> net;
	std::vector<std::vector<ID>> link;
	double m_width, m_height;

public:
	RBNet() {}
	RBNet(double _w, double _h, const std::__1::vector<Point> &vec)
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
