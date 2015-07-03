//
// Created by Kanari on 15/7/3.
//

#include <vector>
#include "RBRouter.h"

using namespace std;

double RBRouter::solve() {

	// Net ordering
	vector<unsigned int> seq = this->find_order(net);
	// Sequential net embedding
	RBSequentialEmbedding embedding(net, seq);
	// ROAR optimization
	for (auto it = embedding.point_iterator_begin();
		 it != embedding.point_iterator_end(); ++it) {
		embedding.roar(it);
	}

	double result = embedding.length();
	return result;
}

pair<double, RoutingPlan *> RBRouter::solve_with_plan() {
	// Unimplemented, leave here
	return pair<double, RoutingPlan *>();
}

vector<unsigned int> RBRouter::find_order(const RBNet &net) {
	unsigned int n = net.n_nets();
	vector<unsigned int> seq(n);
	for (unsigned int i = 0; i < n; ++i)
		seq[i] = i;

	// Implementation here

	return seq;
}

