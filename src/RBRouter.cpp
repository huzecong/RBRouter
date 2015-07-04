//
// Created by Kanari on 15/7/3.
//

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
	// ROAR optimization
	for (auto it = embedding.point_iterator_begin();
		 it != embedding.point_iterator_end(); ++it) {
		embedding.roar(it);
	}

	double result = embedding.length();
	this->plan = embedding.routing_plan();
	return result;
}

vector<unsigned int> RBRouter::find_order(const RBNet &net) {
	unsigned int max = net.n_nets(); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	vector<unsigned int> seq(max);
	for (unsigned int i = 0; i < max; ++i)
		seq[i] = i;
	unsigned int m = net.n_points();
	
	// Implementation here
	typedef unsigned int ID;
	
	// First find all of the components
	std::vector < std::vector<ID> > cpnID;
	std::vector<RBNet> components;
	std::vector<double> lengths;
	bool *is_linked = new bool[m];
	for(int i = 0; i < m; ++i)
		is_linked[i] = 0;

	debug("\tfind_order: init");
	
	for(int i = 0; i < m; ++i) {
		if(!is_linked[i]) {
			std::queue<ID> q;
			std::vector<ID> cpn;
			q.push(i);
			int y = i;
			while(!q.empty()) {
				if(!is_linked[y]) {
					is_linked[y] = 1;
					cpn.push_back(y);
					std::vector<ID> lfi = net.links_from(y);
					for(unsigned ii = 0; ii < lfi.size(); ++ii) {
						q.push(lfi[ii]);
					}
				}
				q.pop();
				y = q.front();
			}
			cpnID.push_back(cpn);
			debug("\t\tfind_order: before subnet " << i);
			RBNet rb = net.subnet(cpn);
			debug("\t\tfind_order: after subnet " << i);
			components.push_back(rb);
			debug("\t\tfind_order: before embedding " << i);
			RBSequentialEmbedding re(rb, cpn);
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
			
			std::vector<unsigned int> _seq;
			_seq = cpnID[i];
			for(int ii = 0; ii < cpnID[j].size(); ++ii)
				_seq.push_back(cpnID[j][ii]);
			RBNet _net = net.subnet(_seq);	
			
			RBSequentialEmbedding aij(_net, _seq);
			Matrix[i][j] = aij.length() - lengths[i] - lengths[j];
			_s += Matrix[i][j];
		}
		s[i] = _s;
	}
	
	std::vector <ID> ret(n);
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
	
	for(int i = 0; i < n; ++i) {
		delete [] s;
		delete [] selected;
		delete [] Matrix[i];
	}
	delete [] Matrix;
	return seq;
}

void RBRouter::plot(string filename) const {
	vector<PSLine> lines;
	for (const vector<ID> &vec : this->plan.path) {
		PSColor color;
		color.r_ = rand() % 256;
		color.g_ = rand() % 256;
		color.b_ = rand() % 256;
		color.alpha_ = 255;
		for (int i = 0; i + 1 < vec.size(); ++i) {
			Point a = net.point[vec[i]];
			Point b = net.point[vec[i + 1]];
			PSLine line;
			line.x1_ = a.x;
			line.y1_ = a.y;
			line.x2_ = b.x;
			line.y2_ = b.y;
			line.width_ = 1.0;
			line.color_ = color;
			lines.push_back(line);
		}
	}

	PSPlot plot;
	plot.setLines(lines);
	plot.draw(filename);
}
