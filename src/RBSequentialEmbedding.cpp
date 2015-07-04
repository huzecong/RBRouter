//
// Created by Kanari on 15/7/3.
//

#include <climits>
#include <functional>
#include <utility>
#include <vector>
#include "RBSequentialEmbedding.h"
#if __RB_DEBUG
#include "HYPlotPS/HYPlotPS.h"
#endif

using namespace std;

ID RBSERegion::id_cnt;
vector<RBSERegion *> RBSERegion::regions;
vector<ID> RBSERegion::free_id;
b2BlockAllocator *RBSERegion::allocator;
b2BlockAllocator *RBSEVertex::allocator;
RBSEPort RBSEPort::PortNotFound(DBL_MAX, DBL_MAX);

inline RBSEPort make_port(double s, double e) {
	return RBSEPort(s, e);
}

inline double angle(RBSERegion *a, RBSERegion *b) {
	return angle(a->belong_vertex->point, b->belong_vertex->point);
}

bool RBSEPort::contains_eq(double angle) const {
	if (this->start_angle >= this->end_angle) {
		if (this->start_angle <= angle && angle <= this->end_angle + 2 * M_PI)
			return true;
		if (this->start_angle - 2 * M_PI <= angle && angle <= this->end_angle)
			return true;
		return false;
	} else return this->start_angle <= angle && angle <= this->end_angle;
}

bool RBSEPort::contains(double angle) const {
	if (this->start_angle >= this->end_angle) {
		if (lt(this->start_angle, angle)
			&& lt(angle, this->end_angle + 2 * M_PI))
			return true;
		if (lt(this->start_angle - 2 * M_PI, angle)
			&& lt(angle, this->end_angle))
			return true;
		return false;
	} else return lt(this->start_angle, angle)
				  && lt(angle, this->end_angle);
}

bool RBSEPort::contains(const RBSEPort &other) const {
	return this->contains(other.start_angle) && this->contains(other.end_angle)
		   && ~((this->start_angle >= this->end_angle)
				^ (other.start_angle == other.end_angle));
}

double RBSEPort::length() const {
	if (this->start_angle >= this->end_angle) {
		return this->end_angle - this->start_angle + 2 * M_PI;
	} else return this->end_angle - this->start_angle;
}

// CCW order
template<typename T>
T find_prev(LinkedList<T> &list, double angle) {
	const auto func = [angle](const T &data) -> double {
		double t_angle = data->angle;
		if (t_angle > angle) t_angle -= M_PI;
//		if (angle - t_angle <= M_PI) {
			return angle - t_angle;
//		} else return DBL_MAX;
	};
	auto *node = list.find_minimum(func);
	if (node == NULL) return T();
	else return node->data;
}

template<typename T>
T find_next(LinkedList<T> &list, double angle) {
	const auto func = [angle](const T &data) -> double {
		double t_angle = data->angle;
		if (t_angle < angle) t_angle += M_PI;
//		if (t_angle - angle <= M_PI) {
			return t_angle - angle;
//		} else return DBL_MAX;
	};
	auto *node = list.find_minimum(func);
	if (node == NULL) return T();
	else return node->data;
}

RBSEPort find_port(LinkedList<RBSEPort> &list, double angle) {
	const auto func = [angle](const RBSEPort &data) -> double {
		return data.contains(angle);
	};
	auto node = list.find_if(func);
	if (node == NULL) return RBSEPort::PortNotFound;
	else return node->data;
}

void insert_net(LinkedList<RBSENet *> &list, RBSENet *net) {
	double angle = net->angle;
	// find_next
	const auto func = [angle](RBSENet * const &data) -> double {
		double t_angle = data->angle;
		if (t_angle < angle) t_angle += M_PI;
		if (t_angle - angle <= M_PI) {
			return t_angle - angle;
		} else return DBL_MAX;
	};
	auto *x = list.find_minimum(func);
	if (x == NULL) list.append(net);
	else list.append(x, net);
}


RBSequentialEmbedding::RBSequentialEmbedding(const RBNet &net,
											 const vector<unsigned int> &seq) {
	debug("Before embedding");
#if __RB_DEBUG
	debug("\tembed sequence:");
	cerr << "\t";
	for (auto x : seq)
		cerr << x << " ";
	cerr << endl;
#endif
	// Initialize
	this->net = net;
	RBSEVertex::allocator = new b2BlockAllocator();
	RBSERegion::allocator = new b2BlockAllocator();
	RBSERegion::init();

	// Create initial open regions and add all pairs of edges
	for (int i = 0; i < net.n_points(); ++i) {
		RBSEVertex *v = new RBSEVertex();
		v->id = i;
		v->point = net.point[i];
		RBSERegion *region = new RBSERegion();
		region->belong_vertex = v;
		region->port.append(make_port(-M_PI, M_PI));
		region->is_open = true;
		for (int j = 0; j < i; ++j) {
			bool valid = true;
			for (int k = 0; k < net.n_points() && valid; ++k)
				if (k != i && k != j) {
					if (colinear(net.point[i], net.point[j], net.point[k]))
						valid = false;
				}
			if (!valid) continue;
			RBSERegion *to = vertex[j]->region_list.front()->data;
			ID e = graph.add_edge(net.n_points() + region->id,
								  net.n_points() + to->id,
								  dist(vertex[j]->point, v->point) - eps);
			region->link.append(make_pair(to, e));
		}
		ID e = graph.add_edge(i, net.n_points() + region->id, eps);
		v->open_region.append(make_pair(region, e));
		v->region_list.append(region);
		this->vertex.push_back(v);
	}

	debug("\tembed: init");

	// Embed nets one by one
	for (int _net = 0; _net < net.n_nets(); ++_net) {
		// Find shortest path of regions
		ID net_id = seq[_net];
		pair<ID, ID> cur_net = net.net[net_id];
		vector<ID> path = graph.shortest_path(cur_net.first, cur_net.second);
		vector<ID> plan_path;
		for (int i = 1; i + 1 < path.size(); ++i)
			path[i] -= net.n_points();
		for (int i = 1; i + 1 < path.size(); ++i)
			plan_path.push_back(RBSERegion::regions[path[i]]
										->belong_vertex->id);
		this->plan.path.push_back(plan_path);

		// Split regions along the path
		// First and last nodes on path are nodes for vertexes, not regions
		for (int i = 1; i + 1 < path.size(); ++i) {
			ID x = path[i];
			RBSERegion *region = RBSERegion::regions[x];
			RBSERegion *nr[2] = { NULL, NULL };

			// Clear original edges
			for (auto *p = region->link.head;
				 p != region->link.tail; p = p->next) {
				graph.remove_edge(p->data.second);
				p->data.first->link.remove(make_pair(region, p->data.second));
			}
			if (region->is_open) {
				auto *x = region->belong_vertex->open_region.
						find_if([region](const pair<RBSERegion *, ID> &data) {
					return data.first == region;
				});
				graph.remove_edge(x->data.second);
				region->belong_vertex->open_region.remove(x);
			}
			region->belong_vertex->region_list.remove(region);

			if (i == 1 || i + 2 == path.size()) {
				// Incident net
				RBSEIncidentNet *net = new RBSEIncidentNet();
				net->counterpart = NULL;
				if (i == 1) net->from_node = path[i + 1];
				else net->from_node = path[i - 1];
				net->net_no = net_id;
				net->angle = angle(region, RBSERegion::regions[net->from_node]);

				// Split the region
				RBSENet *prev = find_prev(region->net, net->angle);
				RBSENet *next = find_next(region->net, net->angle);
				assert(prev == NULL && next == NULL
					   || prev != NULL && next != NULL);
				assert(region->is_open == true);
				if (prev == NULL && next == NULL) {
					// First net inserted
					// Simply change port
					nr[0] = new RBSERegion();
					nr[0]->port.append(make_port(net->angle, net->angle));
					insert_net(nr[0]->net, net);
					nr[0]->is_open = true;
				} else if (prev->type == AttachedNet
						   && next->type == AttachedNet) {
					// Is not adjacent to any incident nets
					// Modify original region
					nr[0] = new RBSERegion(*region);
					RBSEPort port = find_port(nr[0]->port, net->angle);
					nr[0]->port.remove(port);
					nr[0]->port.append(make_port(port.start_angle,
												 net->angle));
					nr[0]->port.append(make_port(net->angle,
												 port.end_angle));
					insert_net(nr[0]->net, net);
					nr[0]->is_open = region->is_open;
				} else {
					if (prev->type == IncidentNet) {
						// Create new region with left one
						nr[0] = new RBSERegion();
						nr[0]->port.append(make_port(prev->angle, net->angle));
						insert_net(nr[0]->net, prev);
						insert_net(nr[0]->net, net);
						nr[0]->is_open = true;
					} else {
						// Modify original region for left
						nr[0] = new RBSERegion(*region);
						RBSEPort port = find_port(nr[0]->port, net->angle);
						nr[0]->port.remove(port);
						nr[0]->port.append(make_port(net->angle,
													 port.end_angle));
					}
					if (next->type == IncidentNet) {
						// Create new region with right one
						nr[1] = new RBSERegion();
						nr[1]->port.append(make_port(net->angle, next->angle));
						insert_net(nr[1]->net, net);
						insert_net(nr[1]->net, next);
						nr[0]->is_open = true;
					} else {
						// Modify original region for right
						nr[1] = new RBSERegion(*region);
						RBSEPort port = find_port(nr[1]->port, net->angle);
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

				// Split the region
				RBSEPort port[2];
				port[0] = find_port(region->port, net[0]->angle);
				port[1] = find_port(region->port, net[1]->angle);
				if (port[0] == port[1]) {
					// Creates an one-port region
					nr[0] = new RBSERegion();
					nr[0]->is_open = false;
					insert_net(nr[0]->net, net[0]);
					insert_net(nr[0]->net, net[1]);
					nr[0]->port.append(make_port(net[0]->angle, net[1]->angle));
					nr[1] = new RBSERegion(*region);
					insert_net(nr[1]->net, net[0]);
					insert_net(nr[1]->net, net[1]);
					nr[1]->port.remove(port[0]);
					nr[1]->port.append(make_port(port[0].start_angle,
												 net[0]->angle));
					nr[1]->port.append(make_port(net[1]->angle,
												 port[0].end_angle));
				} else {
					// "Horizontally" split the region into two
					RBSENet *middle[2];
					middle[0] = region->net.find_if([&](const RBSENet *x) {
						return equal(x->angle, port[0].end_angle);
					})->data;
					middle[1] = region->net.find_if([&](const RBSENet *x) {
						return equal(x->angle, port[1].start_angle);
					})->data;
					if (middle[0] == middle[1]) {
						// Separated by one incident net
						assert(middle[0]->type == IncidentNet);
					} else {
						// Separated by two incident nets
						// or two parts of one attached net
						if (find_next(region->net, middle[0]->angle)
							!= middle[1]) {
							middle[0] = region->net.find_if([&](const RBSENet *x) {
								return equal(x->angle, port[1].end_angle);
							})->data;
							middle[1] = region->net.find_if([&](const RBSENet *x) {
								return equal(x->angle, port[0].start_angle);
							})->data;
							assert(find_next(region->net, middle[0]->angle) \
								   == middle[1]);
							swap(net[0], net[1]);
						}
						assert(middle[0]->type == middle[1]->type);
						if (middle[0]->type == AttachedNet) {
							assert(middle[0]->counterpart == middle[1]);
						}
					}

					nr[0] = new RBSERegion();
					nr[0]->is_open = false;
					insert_net(nr[0]->net, net[0]);
					insert_net(nr[0]->net, net[1]);
					insert_net(nr[0]->net, middle[0]);
					if (middle[0] != middle[1])
						insert_net(nr[0]->net, middle[1]);
					nr[0]->port.append(make_port(net[0]->angle,
												 port[0].end_angle));
					nr[0]->port.append(make_port(port[1].start_angle,
												 net[1]->angle));

					nr[1] = new RBSERegion(*region);
					nr[1]->net.remove(middle[0]);
					if (middle[0] != middle[1])
						nr[1]->net.remove(middle[1]);
					insert_net(nr[1]->net, net[0]);
					insert_net(nr[1]->net, net[1]);
					nr[1]->port.remove(port[0]);
					nr[1]->port.remove(port[1]);
					nr[1]->port.append(make_port(port[0].start_angle,
												 net[0]->angle));
					nr[1]->port.append(make_port(net[1]->angle,
												 port[1].end_angle));
				}
			}

			// Validate and add edges for new regions
			for (int _nr = 0; _nr < 2; ++_nr) {
				RBSERegion *new_region = nr[_nr];
				if (new_region == NULL) break;
				new_region->belong_vertex = region->belong_vertex;
				new_region->link.reset();
				for (auto *p = region->link.head;
					 p != region->link.tail; p = p->next) {
					RBSERegion *node = p->data.first;
					RBSEPort port = find_port(new_region->port,
											  angle(new_region, node));
					ID x = graph.add_edge(net.n_points() + new_region->id,
										  net.n_points() + node->id,
										  dist(new_region->belong_vertex->point,
											   node->belong_vertex->point));
					new_region->link.append(make_pair(node, x));
				}
				region->belong_vertex->region_list.append(new_region);
				if (new_region->is_open) {
					ID e = graph.add_edge(region->belong_vertex->id,
										  net.n_points() + new_region->id, eps);
					region->belong_vertex->
							open_region.append(make_pair(new_region, e));
				}
			}
		}
		for (int i = 1; i + 1 < path.size(); ++i)
			delete RBSERegion::regions[path[i]];

		debug("\t\tembed: finished net " << _net);
	}

	// Compute total length
	this->plan.length = 0.0;
	for (const vector<ID> &vec : this->plan.path) {
		for (int i = 0; i + 1 < vec.size(); ++i)
			this->plan.length += dist(net.point[vec[i]], net.point[vec[i + 1]]);
	}
	debug("\tembed length: " << this->plan.length);
	debug("After embedding");

	static char run_times = 'a' - 1;
	++run_times;
	string filename = "run";
	filename.append(1, run_times);
	filename += ".ps";
	plot_postscript(filename, this->net.point, this->plan);
}

void RBSequentialEmbedding::roar(const Point_Iterator &it) {

}
