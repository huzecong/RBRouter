//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_SEQUENTIALEMBEDDING_H
#define RBROUTER_SEQUENTIALEMBEDDING_H


#include <vector>
#include <utility>
#include "RBBase.h"


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
			  link(vec.size(), std::__1::vector<ID>()) {}
	ID add_point(const Point &p) {
		point.push_back(p);
		link.push_back(std::__1::vector<ID>());
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
	RBNet subnet(std::__1::vector<ID> vec) const;
};

class RBSEAttachedNet {

};

class RBSEIncidentNet {

};

class RBSEPath {

};

class RBSEPoint {
	RBSEPath path;

};

class RBSequentialEmbedding {
	typedef std::vector<RBSEPoint>::iterator Point_Iterator;
	RBNet net;
	std::vector<RBSEPoint> point;

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
