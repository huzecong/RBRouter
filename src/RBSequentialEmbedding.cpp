//
// Created by Kanari on 15/7/3.
//

#include <functional>
#include <utility>
#include <vector>
#include "RBSequentialEmbedding.h"

using namespace std;

ID RBSERegion::id_cnt;
vector<RBSERegion *> RBSERegion::regions;
vector<ID> RBSERegion::free_id;
b2BlockAllocator *RBSERegion::allocator;
b2BlockAllocator *RBSEVertex::allocator;

inline RBSEPort make_port(double s, double e) {
	return RBSEPort(s, e);
}

inline double angle(RBSERegion *a, RBSERegion *b) {
	return angle(a->belong_vertex->point, b->belong_vertex->point);
}

inline bool port_contains_eq(const RBSEPort &port, double angle) {
	if (port.start_angle >= port.end_angle) {
		if (port.start_angle <= angle && angle <= port.end_angle + 2 * M_PI)
			return true;
		if (port.start_angle - 2 * M_PI <= angle && angle <= port.end_angle)
			return true;
		return false;
	} else return port.start_angle <= angle && angle <= port.end_angle;
}

inline bool port_contains(const RBSEPort &port, double angle) {
	if (port.start_angle >= port.end_angle) {
		if (less(port.start_angle, angle)
			&& less(angle, port.end_angle + 2 * M_PI))
			return true;
		if (less(port.start_angle - 2 * M_PI, angle)
			&& less(angle, port.end_angle))
			return true;
		return false;
	} else return less(port.start_angle, angle) && less(angle, port.end_angle);
}

inline bool port_contains(const RBSEPort &port, const RBSEPort &other) {
	return port_contains(port, other.start_angle)
		   && port_contains(port, other.end_angle)
		   && ~((port.start_angle >= port.end_angle)
				^ (other.start_angle == other.end_angle));
}

// CCW order
template<typename T>
T find_prev(LinkedList<T> &list, RBSENet *net) {
	const auto func = [&](const T &data) -> double {
		double angle = data->angle;
		if (angle > net->angle) angle -= M_PI;
		if (net->angle - angle <= M_PI) {
			return net->angle - angle;
		} else return DBL_MAX;
	};
	return list.find_minimum(func)->data;
}

template<typename T>
T find_next(LinkedList<T> &list, RBSENet *net) {
	const auto func = [&](const T &data) -> double {
		double angle = data->angle;
		if (angle < net->angle) angle += M_PI;
		if (angle - net->angle <= M_PI) {
			return angle - net->angle;
		} else return DBL_MAX;
	};
	return list.find_minimum(func)->data;
}

RBSEPort find_port(LinkedList<RBSEPort> &list, RBSENet *net) {
	const auto func = [&](const RBSEPort &data) -> double {
		return port_contains(data, net->angle);
	};
	return list.find_if(func)->data;
}


RBSequentialEmbedding::RBSequentialEmbedding(const RBNet &net,
											 const vector<unsigned int> &seq) {
	// Initialize
	this->net = net;
	RBSEVertex::allocator = new b2BlockAllocator();
	RBSERegion::allocator = new b2BlockAllocator();

	// Create initial open regions and add all pairs of edges
	for (int i = 0; i < net.n_points(); ++i) {
		RBSEVertex *v = new RBSEVertex();
		v->point = net.point[i];
		RBSERegion *region = new RBSERegion();
		region->belong_vertex = v;
		region->port.append(make_port(0.0, 2 * M_PI));
		region->is_open = true;
		for (int j = 0; j < i; ++j) {
			RBSERegion *to = vertex[j]->region_list.front()->data;
			ID e = graph.add_edge(region->id, to->id, dist(vertex[j]->point,
														   v->point));
			region->link.append(make_pair(to, e));
		}
		v->open_region.append(region);
		v->region_list.append(region);
		this->vertex.push_back(v);
	}

	// Embed nets one by one
	for (int _net = 0; _net < net.n_nets(); ++_net) {
		// Find shortest path of regions
		ID net_id = seq[_net];
		pair<ID, ID> cur_net = net.net[net_id];
		vector<ID> path = graph.shortest_path(cur_net.first, cur_net.second);
		plan.path.push_back(path);

		// Split regions along the path
		for (int i = 0; i < path.size(); ++i) {
			ID x = path[i];
			RBSERegion *region = RBSERegion::regions[x];
			RBSERegion *nr[2] = { NULL, NULL };
			region->belong_vertex->region_list.remove(region);
			region->belong_vertex->open_region.remove(region);
			if (i == 0 || i + 1 == path.size()) {
				// Incident net
				RBSEIncidentNet *net = new RBSEIncidentNet();
				net->counterpart = NULL;
				if (i == 0) net->from_node = path[i + 1];
				else net->from_node = path[i - 1];
				net->net_no = net_id;
				net->angle = angle(region, RBSERegion::regions[net->from_node]);

				// Split the region
				RBSENet *prev = find_prev(region->net, net);
				RBSENet *next = find_next(region->net, net);
				assert(prev == NULL && next == NULL
					   || prev != NULL && next != NULL);
				assert(region->is_open == true);
				if (prev == NULL && next == NULL) {
					// First net inserted
					// Simply change port
					nr[0] = new RBSERegion();
					nr[0]->port.append(make_port(net->angle, net->angle));
					nr[0]->net.append(net);
					nr[0]->is_open = true;
				} else if (prev->type == AttachedNet
						   && next->type == AttachedNet) {
					// Is not adjacent to any incident nets
					// Modify original region
					nr[0] = new RBSERegion(*region);
					RBSEPort port = find_port(nr[0]->port, net);
					nr[0]->port.remove(port);
					nr[0]->port.append(make_port(port.start_angle,
												 net->angle));
					nr[0]->port.append(make_port(net->angle,
												 port.end_angle));
					nr[0]->is_open = region->is_open;
				} else {
					if (prev->type == IncidentNet) {
						// Create new region with left one
						nr[0] = new RBSERegion();
						nr[0]->port.append(make_port(prev->angle, net->angle));
						nr[0]->net.append(prev);
						nr[0]->net.append(net);
						nr[0]->is_open = true;
					} else {
						// Modify original region for left
						nr[0] = new RBSERegion(*region);
						RBSEPort port = find_port(nr[0]->port, net);
						nr[0]->port.remove(port);
						nr[0]->port.append(make_port(net->angle,
													 port.end_angle));
					}
					if (next->type == IncidentNet) {
						// Create new region with right one
						nr[1] = new RBSERegion();
						nr[1]->port.append(make_port(net->angle, next->angle));
						nr[1]->net.append(net);
						nr[1]->net.append(next);
						nr[0]->is_open = true;
					} else {
						// Modify original region for right
						nr[1] = new RBSERegion(*region);
						RBSEPort port = find_port(nr[1]->port, net);
						nr[1]->port.remove(port);
						nr[1]->port.append(make_port(port.start_angle,
													 net->angle));
					}
				}

			} else {
				// Attached net
				RBSEAttachedNet *net[2];
				net[0] = new RBSEAttachedNet();
				net[1] = new RBSEAttachedNet();
				net[0]->counterpart = net[1];
				net[1]->counterpart = net[0];
				net[0]->from_node = path[i - 1];
				net[1]->from_node = path[i + 1];
				net[0]->net_no = net[1]->net_no = net_id;
				net[0]->angle = angle(region, RBSERegion::regions[path[i - 1]]);
				net[1]->angle = angle(region, RBSERegion::regions[path[i + 1]]);

				RBSEPort port[2];
				port[0] = find_port(region->port, net[0]);
				port[1] = find_port(region->port, net[1]);
				if (port[0] == port[1]) {

				}
			}

			// Validate and add edges for new regions
			for (int _nr = 0; _nr < 2; ++_nr) {
				RBSERegion *new_region = nr[_nr];
				new_region->belong_vertex = region->belong_vertex;
			}
		}
	}

}

void RBSequentialEmbedding::roar(const Point_Iterator &it) {

}

double RBSequentialEmbedding::length() {
	return 0;
}
