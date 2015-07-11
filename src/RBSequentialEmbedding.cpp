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
using namespace std::rel_ops;

ID RBSERegion::id_cnt;
vector<RBSERegion *> RBSERegion::regions;
vector<ID> RBSERegion::free_id;
b2BlockAllocator *RBSERegion::allocator;
b2BlockAllocator *RBSEVertex::allocator;

inline RBSEPort make_port(RBSENet *l, RBSENet *r) {
	RBSEPort ret = RBSEPort(l->to_vertex->point - l->from_vertex->point,
							r->to_vertex->point - r->from_vertex->point);
//	assert(sign(cross_product(ret.s, ret.e)) >= 0);
	return ret;
}

inline double angle(RBSERegion *a, RBSERegion *b) {
	return angle(a->vertex->point, b->vertex->point);
}

inline double reverse(double x) {
	return x > 0 ? x -= M_PI : x + M_PI;
}

inline Point reverse(const Point &p) {
	return Point(-p.x, -p.y);
}

inline bool turn_left(RBSENet *a, RBSENet *b) {
	assert(a->from_vertex == b->from_vertex);
	return sign(cross_product(a->from_vertex->point,
							  a->to_vertex->point,
							  b->to_vertex->point)) <= 0;
}

enum RBSERegionNetSide {
	NotAdjacent = 0,
	LeftSide = 1,
	RightSide = 2,
	BothSides = 3
};

// Return on which side of the net is the region
inline RBSERegionNetSide on_side_of(RBSERegion *region, RBSENet *net) {
	if (net->type == AttachedNet) {
		if (region->outer_net[0] == net) return RightSide;
		if (region->outer_net[1] == net) return LeftSide;
		if (region->inner_net[0] == net) return LeftSide;
		if (region->inner_net[1] == net) return RightSide;
		return NotAdjacent;
	} else {
		if (region->incident_net[0] == net) {
			if (region->incident_net[1] == net) return BothSides;
			return LeftSide;
		} else {
			if (region->incident_net[1] == net) return RightSide;
			return NotAdjacent;
		}
	}
}

bool same_way(const Point &a, const Point &b) {
	if (sign(cross_product(a, b)) != 0) return false;
	return sign(a.x) == sign(b.x) && sign(a.y) == sign(b.y);
}

bool RBSEPort::contains_eq(const Point &p) const {
	if (this->contains(p)) return true;
	if (same_way(this->s, p)) return true;
	if (same_way(this->e, p)) return true;
	return false;
}

bool RBSEPort::contains(const Point &p) const {
	int cps = sign(cross_product(this->s, p));
	int cpe = sign(cross_product(p, this->e));
	if (cps > 0 && cpe > 0) return true;
	if (sign(cross_product(this->s, this->e)) >= 0) return false;
	if (cpe > 0) {
		cps = sign(cross_product(reverse(this->s), p));
	} else if (cps > 0) {
		cpe = sign(cross_product(p, reverse(this->e)));
	}
	return cps > 0 && cpe > 0;
}

// CCW order
template<typename T>
T find_prev(LinkedList<T> &list, double angle) {
	const auto func = [angle](const T &data) -> double {
		double t_angle = data->angle;
		if (t_angle > angle) t_angle -= 2 * M_PI;
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
		if (t_angle < angle) t_angle += 2 * M_PI;
//		if (t_angle - angle <= M_PI) {
			return t_angle - angle;
//		} else return DBL_MAX;
	};
	auto *node = list.find_minimum(func);
	if (node == NULL) return T();
	else return node->data;
}

// Insert before
inline void insert_net(LinkedList<RBSENet *> &list, RBSENet *x, RBSENet *net) {
	auto *node = list.find(x);
	list.append(node, net);
}

inline void insert_net(LinkedList<RBSENet *> &list, RBSENet *net) {
	list.append(net);
}

inline bool intersect(RBSERegion *a, RBSERegion *b,
					  RBSERegion *s, RBSERegion *t) {
	return intersect_strict(a->vertex->point,
							b->vertex->point,
							s->vertex->point,
							t->vertex->point);
}

vector<ID> RBSequentialEmbedding::insert_open_edge(RBSEVertex *vertex) {
	vector<ID> e;
	for (auto *p = vertex->region_list.head;
		 p != vertex->region_list.tail; p = p->next) {
		if (p->data->incident_net[0] == NULL
			&& (p->data->outer_net[0] != NULL || p->data->inner_net[0] != NULL))
			continue;
		ID x = graph.add_edge(vertex->id, p->data->id + net.n_points(), eps);
		e.push_back(x);
	}
	return e;
}


void RBSequentialEmbedding::preprocess() {
	// Add corner vertex
	ID lu = net.n_points(), ru = lu + 1, rd = ru + 1, ld = rd + 1;
	for (unsigned i = 0; i < net.n_points(); ++i) {
		Point &p = net.point[i];
		if (equal(p.x, 0.0)) {
			if (equal(p.y, 0.0)) ld = i;
			else if (equal(p.y, net.height())) lu = i;
		} else if (equal(p.x, net.width())) {
			if (equal(p.y, 0.0)) rd = i;
			else if (equal(p.y, net.height())) ru = i;
		}
	}
	if (lu >= net.n_points())
		net.add_point(Point(0, net.height()));
	if (ru >= net.n_points())
		net.add_point(Point(net.width(), net.height()));
	if (rd >= net.n_points())
		net.add_point(Point(net.width(), 0));
	if (ld >= net.n_points())
		net.add_point(Point(0, 0));

	// Create initial open regions and add all pairs of edges
	for (int i = 0; i < net.n_points(); ++i) {
		RBSEVertex *v = new RBSEVertex();
		v->id = i;
		v->point = net.point[i];
		RBSERegion *region = new RBSERegion();
		region->vertex = v;
		for (int j = 0; j < i; ++j) {
			bool valid = true;
			for (int k = 0; k < net.n_points() && valid; ++k)
				if (k != i && k != j
					&& between(net.point[k], net.point[i], net.point[j])) {
					if (colinear(net.point[i], net.point[j], net.point[k]))
						valid = false;
				}
			if (!valid) continue;
			RBSERegion *to = vertex[j]->region_list.front()->data;
			ID e = graph.add_edge(net.n_points() + region->id,
								  net.n_points() + to->id,
								  dist(vertex[j]->point, v->point) - eps);
			region->link.append(make_pair(to, e));
			to->link.append(make_pair(region, e));
		}
		v->region_list.append(region);
		this->vertex.push_back(v);
	}

	// Process vertex on border
	vector<ID> vec[4];
	for (unsigned i = 0; i < net.n_points(); ++i) {
		Point &p = net.point[i];
		if (equal(p.x, 0.0)) vec[0].push_back(i);
		if (equal(p.x, net.width())) vec[1].push_back(i);
		if (equal(p.y, 0.0)) vec[2].push_back(i);
		if (equal(p.y, net.height())) vec[3].push_back(i);
	}
	sort(vec[0].begin(), vec[0].end(), [&](int x, int y) {
		return net.point[x].y < net.point[y].y;
	});
	sort(vec[1].begin(), vec[1].end(), [&](int x, int y) {
		return net.point[x].y < net.point[y].y;
	});
	sort(vec[2].begin(), vec[2].end(), [&](int x, int y) {
		return net.point[x].x < net.point[y].x;
	});
	sort(vec[3].begin(), vec[3].end(), [&](int x, int y) {
		return net.point[x].x < net.point[y].x;
	});
	for (int x = 0; x < 4; ++x) {
		for (int i = 0; i + 1 < vec[x].size(); ++i) {
			RBSERegion *a = vertex[vec[x][i]]->region_list.head->data;
			RBSERegion *b = vertex[vec[x][i + 1]]->region_list.head->data;
			RBSEIncidentNet *net[2];
			net[0] = new RBSEIncidentNet();
			net[0]->net_no = -1;
			net[0]->from_vertex = a->vertex;
			net[0]->to_vertex = b->vertex;
			net[1] = new RBSEIncidentNet();
			net[1]->net_no = -1;
			net[1]->from_vertex = b->vertex;
			net[1]->to_vertex = a->vertex;
			net[0]->opposite = net[1];
			net[1]->opposite = net[0];
			if (a->incident_net[0] != NULL) {
				a->incident_net[1] = net[0];
			} else a->incident_net[0] = net[0];
			insert_net(a->vertex->net_list, net[0]);
			if (b->incident_net[0] != NULL) {
				b->incident_net[1] = net[1];
			} else b->incident_net[0] = net[1];
			insert_net(b->vertex->net_list, net[1]);
		}
	}
	for (int i = 0; i < net.n_points(); ++i) {
		RBSERegion *region = vertex[i]->region_list.head->data;
		Point &p = net.point[i];
		if (equal(p.x, 0.0) || equal(p.x, net.width())) {
			assert(region->incident_net[0] != NULL
				   && region->incident_net[1] != NULL);
			if (equal(p.y, 0.0) || equal(p.y, net.height())) {
				// Corner point
				if (!turn_left(region->incident_net[1], region->incident_net[0]))
					swap(region->incident_net[0], region->incident_net[1]);
			} else {
				bool left_side = equal(p.x, 0.0);
				bool left_port = region->incident_net[0]->to_vertex->point.y
								 < region->incident_net[1]->to_vertex->point.y;
				if (left_side ^ left_port)
					swap(region->incident_net[0], region->incident_net[1]);
			}
		} else {
			if (equal(p.y, 0.0) || equal(p.y, net.height())) {
				assert(region->incident_net[0] != NULL
					   && region->incident_net[1] != NULL);
				bool down_side = equal(p.y, 0.0);
				bool up_port = region->incident_net[0]->to_vertex->point.x
							   < region->incident_net[1]->to_vertex->point.x;
				if (!(down_side ^ up_port))
					swap(region->incident_net[0], region->incident_net[1]);
			}
		}
	}
}

RBSequentialEmbedding::RBSequentialEmbedding(const RBNet &_net,
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
	net = _net;
	RBSEVertex::allocator = new b2BlockAllocator();
	RBSERegion::allocator = new b2BlockAllocator();
	RBSERegion::init();
	this->plan.width = net.width();
	this->plan.height = net.height();

	this->preprocess();

	debug("\tembed: init");

	// Embed nets one by one
	for (int _net = 0; _net < net.n_nets(); ++_net) {
		debug("\t\tembed: start net " << _net << " (" << seq[_net] << ")");
		// Find shortest path of regions
		ID net_id = seq[_net];
		pair<ID, ID> cur_net = net.net[net_id];
		// Add edges from vertex to open regions
		vector<ID> open_ID[2];
		open_ID[0] = insert_open_edge(vertex[cur_net.first]);
		open_ID[1] = insert_open_edge(vertex[cur_net.second]);
		vector<ID> path = graph.shortest_path(cur_net.first, cur_net.second);
		for (int i = 0; i < 2; ++i)
			for (ID x : open_ID[i])
				graph.remove_edge(x);
		if (path.size() == 0) {
			// No path found, embedding failed
			this->plan.length = INFI;
			return ;
		}
		vector<ID> plan_path;
		for (int i = 1; i + 1 < path.size(); ++i)
			path[i] -= net.n_points();
		for (int i = 1; i + 1 < path.size(); ++i)
			plan_path.push_back(RBSERegion::regions[path[i]]->vertex->id);
		this->plan.path.push_back(plan_path);

		// Remove intersecting edges
		for (auto *region : RBSERegion::regions) {
			if (region == NULL) continue;
			vector<decltype(region->link.head)> marked;
			bool valid = true;
			for (int i = 1; i + 1 < path.size() && valid; ++i) {
				if (region == RBSERegion::regions[path[i]])
					valid = false;
			}
			if (!valid) continue;
			for (auto *x = region->link.head;
				 x != region->link.tail; x = x->next) {
				valid = true;
				for (int i = 1; i + 2 < path.size() && valid; ++i) {
					if (intersect(region, x->data.first,
								  RBSERegion::regions[path[i]],
								  RBSERegion::regions[path[i + 1]]))
						valid = false;
				}
				if (!valid) marked.push_back(x);
			}
			for (auto *x : marked) {
				graph.remove_edge(x->data.second);
				x->data.first->link.remove(make_pair(region, x->data.second));
				region->link.remove(x);
			}
		}


		// Split regions along the path
		// First and last nodes on path are nodes for vertices, not regions
		vector<RBSENet *> last_nets(path.size(), NULL);
		vector<RBSERegion *> new_regions;
		RBSERegion *last_nr[2] = {NULL, NULL};
		vector<int> order(path.size() - 2);
		for (int i = 2; i + 2 < path.size(); ++i)
			order[i - 2] = i;
		order[path.size() - 4] = 1;
		order[path.size() - 3] = path.size() - 2;
		for (int _i = 0; _i < order.size(); ++_i) {
			int i = order[_i];
			ID x = path[i];
			RBSERegion *region = RBSERegion::regions[x];
			RBSERegion *nr[2] = {NULL, NULL};

			// Clear original edges
			for (auto *p = region->link.head; p != region->link.tail; p = p->next) {
				graph.remove_edge(p->data.second);
				p->data.first->link.remove(make_pair(region, p->data.second));
			}
			region->vertex->region_list.remove(region);

			RBSENet *last_net;
			if (i == 1 || i + 2 == path.size()) {
				RBSERegion *other_region;
				if (i == 1) {
					other_region = RBSERegion::regions[path[i + 1]];
					last_net = last_nets[1];
				} else {
					other_region = RBSERegion::regions[path[i - 1]];
					last_net = last_nets[i - 1];
				}
				// Incident net
				RBSEIncidentNet *net = new RBSEIncidentNet();
				net->counterpart = NULL;
				net->net_no = net_id;
				net->from_vertex = region->vertex;
				net->to_vertex = other_region->vertex;
				if (last_net != NULL) {
					net->opposite = last_net;
					last_net->opposite = net;
				}
				last_nets[i] = net;
				if (region->vertex->net_list.size() == 0) {
					// First net inserted
					nr[0] = new RBSERegion(*region);
					nr[0]->incident_net[0] = nr[0]->incident_net[1] = net;
					insert_net(region->vertex->net_list, net);
				} else if (region->outer_net[0] == NULL) {
					// Only incident nets
					nr[0] = new RBSERegion(*region);
					nr[0]->incident_net[0] = net;
					nr[1] = new RBSERegion(*region);
					nr[1]->incident_net[1] = net;
					insert_net(region->vertex->net_list,
							   region->incident_net[1], net);
				} else {
					RBSERegionNetSide res = NotAdjacent;
					Point p = other_region->vertex->point - region->vertex->point;
					RBSEPort port[2];
					port[0] = make_port(region->incident_net[0], region->outer_net[0]);
					port[1] = make_port(region->outer_net[1], region->incident_net[1]);
					if (region->incident_net[0] == region->incident_net[1]
						&& same_way(port[0].s, p)) {
						if (port[0].contains(p)) res = LeftSide;
						else if (port[1].contains(p)) res = RightSide;
						else {
							auto &list = other_region->vertex->net_list;
							auto *x = list.find(region->incident_net[0]->opposite);
							decltype(x) next;
							if (x->next == list.tail) next = list.head;
							else next = x->next;
							if (next->data == last_net) res = LeftSide;
							else res = RightSide;
						}
					} else {
						if (port[0].contains_eq(p)) res = LeftSide;
						else if (port[1].contains_eq(p)) res = RightSide;
					}

					if (res == LeftSide) {
						nr[0] = new RBSERegion(*region);
						nr[0]->outer_net[0] = nr[0]->outer_net[1] = NULL;
						nr[0]->incident_net[1] = net;
						nr[1] = new RBSERegion(*region);
						nr[1]->incident_net[0] = net;
						insert_net(region->vertex->net_list,
								   region->outer_net[0], net);
					} else {
						nr[0] = new RBSERegion(*region);
						nr[0]->outer_net[0] = nr[0]->outer_net[1] = NULL;
						nr[0]->incident_net[0] = net;
						nr[1] = new RBSERegion(*region);
						nr[1]->incident_net[1] = net;
						insert_net(region->vertex->net_list,
								   region->incident_net[1], net);
					}
				}
			} else {
				// Attached net
				last_net = last_nets[i - 1];
				RBSEAttachedNet *net[2];
				net[0] = new RBSEAttachedNet();
				net[1] = new RBSEAttachedNet();
				net[0]->counterpart = net[1];
				net[1]->counterpart = net[0];
				net[0]->net_no = net[1]->net_no = net_id;
				net[0]->from_vertex = net[1]->from_vertex = region->vertex;
				if (last_net == NULL) {
					assert(i == 2);
					last_nets[1] = net[0];
				} else {
					net[0]->opposite = last_net;
					last_net->opposite = net[0];
				}
				last_nets[i] = net[1];
				RBSERegion *other_region[2];
				other_region[0] = RBSERegion::regions[path[i - 1]];
				other_region[1] = RBSERegion::regions[path[i + 1]];
				net[0]->to_vertex = other_region[0]->vertex;
				net[1]->to_vertex = other_region[1]->vertex;
				int res = sign(cross_product(region->vertex->point,
											 other_region[0]->vertex->point,
											 other_region[1]->vertex->point));
				if (res > 0) {
					swap(other_region[0], other_region[1]);
					swap(net[0], net[1]);
				} else if (res == 0
						   && other_region[0]->vertex != other_region[1]->vertex) {
					// not implemented
					assert(false);
				}
				nr[0] = new RBSERegion(*region);
				nr[1] = new RBSERegion(*region);
				nr[0]->incident_net[0] = nr[0]->incident_net[1] = NULL;
				nr[0]->inner_net[0] = net[0];
				nr[0]->inner_net[1] = net[1];
				nr[1]->outer_net[0] = net[0];
				nr[1]->outer_net[1] = net[1];
				if (region->outer_net[0] == NULL) {
					insert_net(region->vertex->attached_list, net[0]);
					if (region->inner_net[0] != NULL) {
						insert_net(region->vertex->net_list, region->inner_net[1], net[1]);
					} else {
						insert_net(region->vertex->net_list, region->incident_net[1], net[1]);
					}
					insert_net(region->vertex->net_list, net[1], net[0]);
				} else {
					insert_net(region->vertex->attached_list, region->outer_net[0], net[0]);
					if (region->inner_net[0] != NULL) {
						insert_net(region->vertex->net_list, region->inner_net[1], net[1]);
					} else {
						insert_net(region->vertex->net_list, region->incident_net[1], net[1]);
					}
					insert_net(region->vertex->net_list, region->outer_net[0], net[0]);
				}
			}

			if (nr[0] != NULL) region->vertex->region_list.append(nr[0]);
			if (nr[1] != NULL) region->vertex->region_list.append(nr[1]);

			for (int _nr = 0; _nr < 2; ++_nr) {
				RBSERegion *new_region = nr[_nr];
				if (new_region == NULL) break;
				new_region->link.reset();
				for (auto *p = region->link.head;
					 p != region->link.tail; p = p->next) {
					RBSERegion *node = p->data.first;
					double w = dist(region->vertex->point, node->vertex->point);
					bool valid = false;
					if (i == 2 && node->id == path[1]) valid = true;
					if (i == 1 && path.size() == 4 && node->id == path[2]) valid = true;
					if (i > 1 && i + 2 < path.size() && node->id == path[i + 1]) valid = true;
					if (new_region->outer_net[0] == NULL && new_region->inner_net[0] == NULL) {
						if (new_region->incident_net[0] == NULL) valid = true;
						if (new_region->incident_net[0] == new_region->incident_net[1])
							valid = true;
					}
					if (!valid) {
						Point p = node->vertex->point - region->vertex->point;
						RBSEPort port[2];
						if (new_region->outer_net[0] == NULL) {
							if (new_region->inner_net[0] == NULL) {
								port[0] = make_port(new_region->incident_net[0],
													new_region->incident_net[1]);
							} else {
								port[0] = make_port(new_region->inner_net[0],
													new_region->inner_net[1]);
							}
							if (port[0].contains(p))
								valid = true;
						} else {
							if (new_region->inner_net[0] == NULL) {
								port[0] = make_port(new_region->incident_net[0],
													new_region->outer_net[0]);
								port[1] = make_port(new_region->outer_net[1],
													new_region->incident_net[1]);
							} else {
								port[0] = make_port(new_region->inner_net[0],
													new_region->outer_net[0]);
								port[1] = make_port(new_region->outer_net[1],
													new_region->inner_net[1]);
							}
							if (port[0].contains(p) || port[1].contains(p))
								valid = true;
						}
					}
					if (!valid) {
						auto on_same_side = [&](RBSERegion *cur_region,
												RBSERegion *other_region,
												RBSENet *this_net) {
							if (this_net == NULL || this_net->opposite == NULL)
								return false;
							auto side_other = on_side_of(other_region, this_net->opposite);
							auto side_this = on_side_of(cur_region, this_net);
							if (side_other == NotAdjacent || side_this == NotAdjacent)
								return false;
							if (side_other == BothSides || side_this == BothSides)
								return true;
							if (side_other != side_this) return true;
							return false;
						};
						for (int i = 0; i < 2; ++i) {
							if (on_same_side(new_region, node, new_region->incident_net[i]))
								valid = true;
							if (on_same_side(new_region, node, new_region->outer_net[i]))
								valid = true;
							if (on_same_side(new_region, node, new_region->inner_net[i]))
								valid = true;
						}
					}
					if (valid) {
						ID x = graph.add_edge(
								net.n_points() + new_region->id,
								net.n_points() + node->id, w);
						new_region->link.append(make_pair(node, x));
						node->link.append(make_pair(new_region, x));
					}
				}
			}

			last_nr[0] = nr[0];
			last_nr[1] = nr[1];
		}

		// Validate and add edges for new regions

		for (int i = 1; i + 1 < path.size(); ++i)
			delete RBSERegion::regions[path[i]];

		debug("\t\tembed: finished net " << _net);
#if __RB_DEBUG
		for (auto *region : RBSERegion::regions) {
			if (region != NULL) {
				debug("\t\t\tcorrespondence " << region->id << ": vertex " <<\
					  region->vertex->id);
			}
		}
#endif
	}

	// Compute total length
	this->plan.length = 0.0;
	for (const vector<ID> &vec : this->plan.path) {
		for (int i = 0; i + 1 < vec.size(); ++i)
			this->plan.length += dist(net.point[vec[i]], net.point[vec[i + 1]]);
	}
	debug("\tembed length: " << this->plan.length);
	debug("After embedding");
/*
	static char run_times = 'a' - 1;
	++run_times;
	string filename = "run";
	filename.append(1, run_times);
	filename += ".ps";
	plot_postscript(filename, this->net.point, this->plan);
*/
}

void RBSequentialEmbedding::roar(ID x) {
	// Not implemented
}


RBSequentialEmbedding::~RBSequentialEmbedding() {

}