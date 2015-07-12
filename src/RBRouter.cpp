//
// Created by Kanari on 15/7/3.
//

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <string>
#include <vector>
#include "RBRouter.h"
#include "HYPlotPS/HYPlotPS.h"

using namespace std;

double RBRouter::solve() {

	this->epsilon_shift();

	// Net ordering
	vector<unsigned int> _seq = this->find_order(net);
	debug("Finished ordering");
#if __RB_DEBUG
	debug("Sequence before PEO:");
	print(_seq);
#endif
	// Planarity enforcement operator
	vector<unsigned int> seq = this->PEO(net, _seq);
#if __RB_DEBUG
	debug("Sequence after PEO:");
	print(seq);
#endif
	// Sequential net embedding
	RBSequentialEmbedding embedding(net, seq);
	debug("Finished embedding");

#if __RB_DEBUG
	debug("Sequence:");
	print(seq);
#endif

	double result = embedding.length();
	this->plan = embedding.routing_plan();
	return result;
}

vector<ID> sequence(int n) {
	vector<ID> ret(n, 0);
	for (int i = 0; i < n; ++i)
		ret[i] = i;
	return ret;
}

vector<ID> RBRouter::find_order(const RBNet &net) {
	unsigned int max = net.n_nets();
	vector<ID> seq = sequence(max);
	unsigned int m = net.n_points();
	
	// Implementation here
	typedef unsigned int ID;
	
	// First find all of the components
	vector<vector<ID>> cpnID;
	vector<RBNet> components;
	vector<double> lengths;
	bool *is_linked = new bool[m];
	for(int i = 0; i < m; ++i)
		is_linked[i] = 0;

	debug("\tfind_order: init m=" << m);
	
	for(int i = 0; i < m; ++i) {
		debug("\t\tfind_order: before checking " << i);
		if(!is_linked[i]) {
			debug("\t\tfind_order: before processing " << i);
			std::queue<ID> q;
			std::vector<ID> nodes, edges;
			q.push(i);
			is_linked[i] = true;
			while(!q.empty()) {
				int y = q.front();
				q.pop();
				nodes.push_back(y);
				const vector<pair<ID, ID>> &lfi = net.links_from(y);
				for(unsigned ii = 0; ii < lfi.size(); ++ii) {
					int node = lfi[ii].first;
					if (node < y) edges.push_back(lfi[ii].second);
					if (is_linked[node]) continue;
					is_linked[node] = true;
					q.push(lfi[ii].first);
					debug("\t\t\t" << y << "--" << lfi[ii].first << ": " << lfi[ii].second);
				}

			}
			cpnID.push_back(edges);
			RBNet rb = net.subnet(nodes);
			components.push_back(rb);
			debug("\t\tfind_order: before embedding " << i);
			debug("\t\tfind_order: |edges|=" << edges.size() << " |nodes|=" << nodes.size());
			vector<unsigned int> seq = this->PEO(rb, sequence(rb.n_nets()));
			RBSequentialEmbedding re(rb, seq);
			debug("\t\tfind_order: after embedding " << i);
			lengths.push_back(re.length());
		}
	}
	delete [] is_linked;

	debug("\tfind_order: components");
	
	// Generate a matrix
	unsigned int n = components.size();
	double** Matrix;
	Matrix = new double*[n];
	bool* selected = new bool[n];
	double* s = new double[n];
	
	for(int i = 0; i < n; ++i) {
		Matrix[i] = new double[n];
		selected[i] = 0;
		double _s = 0;
		for (int j = 0; j < n; ++j) {
			if (i == j) {
				Matrix[i][j] = 0;
				continue;
			}

			RBNet _net = components[i];
			_net.combine(components[j]);

			vector<unsigned int> seq = this->PEO(_net, sequence(_net.n_nets()));
			
			RBSequentialEmbedding aij(_net, seq);
			Matrix[i][j] = aij.length() - lengths[i] - lengths[j];
			_s += Matrix[i][j];
		}
		s[i] = _s;
	}
	
	vector<ID> ret(n);
	for(int i = 0; i < n; ++i) {
		int k = n - i - 1;
		double min = -1;
		unsigned int b = 0;
		for(unsigned int j = 0; j < n; ++j) {
			if(selected[j])
				continue;
			if(min == -1 || min > s[j]) {
				min = s[j];
				b = j;
			}
		}
		selected[b] = 1;
		for(int j = 0; j < n; ++j) {
			s[j] -= Matrix[j][b];
		}
		ret[k] = b;
	}
	reverse(ret.begin(), ret.end());
	
	int counter = 0;
	for(int i = 0; i < ret.size(); ++i) {
		for(int j = 0; j < cpnID[ret[i]].size(); ++j) {
			seq[counter] = cpnID[ret[i]][j];
			++counter;
		}
	}

	debug("\tfind_order: matrix");
#if __RB_DEBUG
	for (int i = 0; i < n; ++i) {
		cerr << "\t\t";
		for (int j = 0; j < n; ++j)
			cerr << Matrix[i][j] << "\t";
		cerr << endl;
	}
#endif

	delete [] s;
	delete [] selected;
	for(int i = 0; i < n; ++i)
		delete [] Matrix[i];
	delete [] Matrix;
	return seq;
}

vector<unsigned int> RBRouter::PEO(const RBNet &net, const vector<unsigned int> &seq) {
	vector<unsigned int> closed, open;
	RBUnionFind ds(net.n_points() + 4);	// 4 Extra points for edges
	ID up = net.n_points(), down = up + 1;
	ID left = down + 1, right = left + 1;
	for (int i = 0; i < net.n_points(); ++i) {
		const Point &p = net.point[i];
		if (equal(p.x, 0.0)) ds.merge(i, left);
		if (equal(p.x, net.width())) ds.merge(i, right);
		if (equal(p.y, 0.0)) ds.merge(i, down);
		if (equal(p.y, net.height())) ds.merge(i, up);
	}
	for (unsigned int x : seq) {
		ID a = net.net[x].first, b = net.net[x].second;
		bool close = false;
		for (int i = 0; i < 2; ++i) {
			if (ds.connected(a, up) && ds.connected(b, down)) close = true;
			if (ds.connected(a, up) && ds.connected(b, left)) close = true;
			if (ds.connected(a, up) && ds.connected(b, right)) close = true;
			if (ds.connected(a, down) && ds.connected(b, left)) close = true;
			if (ds.connected(a, down) && ds.connected(b, right)) close = true;
			if (ds.connected(a, left) && ds.connected(b, right)) close = true;
			swap(a, b);
		}
		if (close) closed.push_back(x);
		else {
			open.push_back(x);
			ds.merge(a, b);
		}
	}
	open.insert(open.end(), closed.begin(), closed.end());
	return open;
}

void RBRouter::plot(string filename) const {
	plot_postscript(filename, this->net.point, this->plan);
}

void RBRouter::epsilon_shift() {
	bool ready = true;
	int adjust_rounds = 0;
	do {
		ready = true;
		double eps_shift = min(net.width(), net.height());
		for (int i = 0; i < net.n_points(); ++i)
			for (int j = i + 1; j < net.n_points(); ++j)
				eps_shift = min(eps_shift, dist(net.point[i], net.point[j]));
		eps_shift /= 100.0;
		for (int i = 0; i < net.n_points(); ++i)
			for (int j = 0; j < net.n_points(); ++j)
				for (int k = 0; k < net.n_points(); ++k) {
					if (i == j || i == k || j == k) continue;
					if (!between(net.point[k], net.point[i], net.point[j]))
						continue;
					if (!colinear(net.point[i], net.point[j], net.point[k]))
						continue;
					Point &p = net.point[k];
					if (equal(p.x, 0.0) || equal(p.x, net.width())
						|| equal(p.y, 0.0) || equal(p.y, net.width()))
						continue;
					ready = false;
					p.x += eps_shift * (rand() % 3 - 1);
					p.y += eps_shift * (rand() % 3 - 1);
				}
		++adjust_rounds;
	} while (!ready);
	debug("eps_shift: " << adjust_rounds << " rounds");
}
