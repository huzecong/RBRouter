//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_RBROUTER_H
#define RBROUTER_RBROUTER_H


#include <string>
#include <utility>
#include <vector>
#include "RBBase.h"
#include "RBSequentialEmbedding.h"

class RBRouter {
public:
	RBNet net;

	std::vector<unsigned int> find_order(const RBNet &);
	std::vector<unsigned int> PEO(const RBNet &, const std::vector<unsigned int> &);

	RBRoutingPlan plan;

public:
	RBRouter() : net() {}
	// Initialize with vector of points, IDs are number from 0
	RBRouter(double w, double h, const std::vector<Point> &vec)
			: net(w, h, vec) {}
	// Returns internal ID of inserted point
	ID insert_point(const Point &p) {
		return net.add_point(p);
	}
	// Insert a net with points specified by their internal IDs
	void insert_net(ID id1, ID id2) {
		net.add_net(id1, id2);
	}
	double &width() { return net.width(); }
	const double &width() const { return net.width(); }
	double &height() { return net.height(); }
	const double &height() const { return net.height(); }

	double solve();
	const RBRoutingPlan &routing_plan() const {
		return plan;
	}

	void plot(std::string) const;
};


#endif //RBROUTER_RBROUTER_H
