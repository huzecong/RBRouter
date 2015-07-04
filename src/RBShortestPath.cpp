//
// Created by Kanari on 15/7/3.
//

#include "RBShortestPath.h"
#include <cstdio>
#include <utility>

using namespace std;

RBShortestPath::RBShortestPath(): n_vertex(0), n_edges(0) {}

ID RBShortestPath::add_edge(ID a, ID b, double w) {
	debug("\t\t\t\tgraph: add edge " << a << " to " << b << " : " << w);
	Edge e = make_pair(n_edges, w);
	Terminal t = make_pair(a, b);
	n_vertex = max(n_vertex, a > b? a + 1: b + 1);
	graph.append(e);
	terminals.append(t);
	++n_edges;
	return n_edges - 1;
}

void RBShortestPath::remove_edge(ID x) {
	LinkedList<Terminal>::ListNode* t = terminals.head;
	for(LinkedList<Edge>::ListNode *p = graph.head;
		p != graph.tail; p = p->next, t = t->next) {
		if((p->data).first == x) {
			debug("\t\t\t\tgraph: remove edge " << t->data.first << " to " << t->data.second << " : " << p->data.second);
			graph.remove(p);
			terminals.remove(t);
		}
	}
}

std::vector<ID> RBShortestPath::shortest_path(ID a, ID b) {
	debug("\t\t\tgraph: shortest path from " << a << " to " << b);

	typedef LinkedList<Terminal>::ListNode TNode;
	typedef LinkedList<Edge>::ListNode ENode;
	std::vector<ID> ret;
	if(a == b) {
		ret.push_back(a);
		return ret;
	}
	double** dist = new double*[n_vertex];
	bool* used = new bool[n_vertex];
	double* length = new double[n_vertex];
	ID* father = new ID[n_vertex];
	
	for(int i = 0; i < n_vertex; ++i) {
		used[i] = 0;
		length[i] = (i == a ? 0 : -1);
		father[i] = (i == a ? a : -1);
		dist[i] = new double[n_vertex];
		for(int j = 0; j < n_vertex; ++j) {
			dist[i][j] = (i == j? 0: -1);
		}
	}

	ENode *e = graph.head;
	for(TNode *t = terminals.head;
		t != terminals.tail; t = t->next, e = e->next) {
		int _a = (t->data).first;
		int _b = (t->data).second;
		double _v = (e->data).second;
		dist[_a][_b] = _v;
		dist[_b][_a] = _v;
	}
	
	int pos = a;
	bool judge = 0;
	for(int i = 0; i < n_vertex; ++i) {
		double w = -1;
		pos = -1;
		for(int j = 0; j < n_vertex; ++j) {
			if(used[j]) continue;
			if(length[j] != -1 && (w == -1 || w > length[j])) {
				pos = j;
				w = length[j];
			}
		}
		if (pos == -1) break;
		used[pos] = 1;
		for(int j = 0; j < n_vertex; ++j) {
			if(j == a || dist[pos][j] == -1 || used[j]) continue;
			if((length[pos] + dist[pos][j] < length[j]) || (length[j] == -1)) {
				length[j] = length[pos] + dist[pos][j];
				father[j] = pos;
			}
		}
	}
	if(length[b] == -1) {
		debug("disconnected!");
		return ret;
	}
	for (int i = 0; i < n_vertex; ++i)
		debug("d[" << i << "] = " << length[i]);

	ID ans = b;
	while(1) {
		ret.insert(ret.begin(), ans);
		if(ans == a) break;
		ans = father[ans];
	}
	
	
	for(int i = 0; i < n_vertex; ++i) {
		delete [] dist[i];
	}
	delete [] dist;
	delete [] length;
	delete [] used;
	delete [] father;
	return ret;
	
}
