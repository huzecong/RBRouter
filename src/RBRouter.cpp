//
// Created by Kanari on 15/7/3.
//

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

	// Net ordering
	vector<unsigned int> seq = this->find_order(net);
	debug("Finished ordering");
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
			std::vector<ID> cpn;
			q.push(i);
			is_linked[i] = true;
			while(!q.empty()) {
				int y = q.front();
				q.pop();
				const vector<pair<ID, ID>> &lfi = net.links_from(y);
				for(unsigned ii = 0; ii < lfi.size(); ++ii) {
					int node = lfi[ii].first;
					if (is_linked[node]) continue;
					is_linked[node] = true;
					q.push(lfi[ii].first);
					cpn.push_back(lfi[ii].second);
				}
			}
			cpnID.push_back(cpn);
			RBNet rb = net.subnet(cpn);
			components.push_back(rb);
			debug("\t\tfind_order: before embedding " << i);
			RBSequentialEmbedding re(rb, sequence(rb.n_nets()));
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
		for(int j = 0; j < n; ++j) {
			if(i == j) {
				Matrix[i][j] = 0;
				continue;
			}

			RBNet _net = components[i];
			_net.combine(components[j]);
			
			RBSequentialEmbedding aij(_net, sequence(_net.n_nets()));
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
	
	int counter = 0;
	for(int i = 0; i < ret.size(); ++i) {
		for(int j = 0; j < cpnID[ret[i]].size(); ++j) {
			seq[counter] = cpnID[ret[i]][j];
			++counter;
		}
	}

	debug("\tfind_order: matrix");

	delete [] s;
	delete [] selected;
	for(int i = 0; i < n; ++i)
		delete [] Matrix[i];
	delete [] Matrix;
	return seq;
}

void RBRouter::plot(string filename) const {
	plot_postscript(filename, this->net.point, this->plan);
}
