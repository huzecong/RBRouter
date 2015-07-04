//
// Created by Kanari on 15/7/3.
//

#include "RBShortestPath.h"
#include <cstdio>
#include <utility>

using namespace std;

RBShortestPath::RBShortestPath(): n_vertex(0), n_edges(0) {}



ID RBShortestPath::add_edge(ID a, ID b, double w) {
	Edge e = make_pair(n_edges, w);
	Terminal t = make_pair(a, b);
	n_vertex = (a > b? a + 1: b + 1);
	if(n_edges == 0) {
		graph.append(graph.head, e);
		terminals.append(terminals.head, t);
	}
	else {
		graph.append(graph.tail, e);
		terminals.append(terminals.tail, t);
	}
	++n_edges;
	return n_edges - 1;
}

void RBShortestPath::remove_edge(ID x) {
	LinkedList<Edge>::ListNode* e_head = graph.head;
	LinkedList<Terminal>::ListNode* t = terminals.head;
	for(LinkedList<Edge>::ListNode *p = graph.head; p != graph.tail; ++p, ++t) {
		if((p->data).first == x) {
			graph.remove(p);
			terminals.remove(t);
		}
	}
}

std::vector<ID> RBShortestPath::shortest_path(ID a, ID b) {
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
	father[a] = a;
	
	for(int i = 0; i < n_vertex; ++i) {
		used[i] = 0;
		length[i] = (i == a? 0: -1);
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
		for(int j = 0; j < n_vertex; ++j) {
			if(used[j]) continue;
			if(w == -1 || (w > length[j] && length[j] != -1)) {
				pos = j;
				w = length[j];
			}
		}
		used[pos] = 1;
		if(pos == b){
			judge = 1;
			break;
		}
		if(!judge) {
			puts("disconnected!");
			return ret;
		}
		for(int j = 0; j < n_vertex; ++j) {
			if(j == a || dist[pos][j] == -1 || used[j]) continue;
			if((length[pos] + dist[pos][j] < length[j]) || (length[j] == -1)) {
				length[j] = length[pos] + dist[pos][j];
				printf("length[%d] = %f\n", pos, length[pos]);
				printf("length[%d] = %f\n", j, length[j]);
				father[j] = pos;
			}
		}
	}
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
