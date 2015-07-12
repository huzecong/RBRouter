//
// Created by Kanari on 15/7/12.
//

#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>
#include "RBBase.h"

using namespace std;

vector<Point> random_points(unsigned int n, double cord = 10.0) {
	mt19937 mt{random_device{}()};
	uniform_real_distribution<double> distribution(0, cord);
	auto rand = [&]() { return distribution(mt); };
	vector<Point> ret;
	for (unsigned int i = 0; i < n; ++i) {
		bool found = false;
		Point p;
		while (!found) {
			p = Point(rand(), rand());
			bool coincide = false;
			for (const Point &q : ret)
				if (p == q) coincide = true;
			if (!coincide) found = true;
		}
		ret.push_back(p);
	}
	return ret;
}

vector<pair<ID, ID>> random_edges(unsigned int n, unsigned int m) {
	mt19937 mt{random_device{}()};
	uniform_int_distribution<int> distribution(0, n - 1);
	auto rand = [&]() { return distribution(mt); };
	set<pair<ID, ID>> h;
	vector<pair<ID, ID>> ret;
	for (int i = 0; i < m; ++i) {
		int a, b;
		do {
			a = rand(), b = rand();
			if (a > b) swap(a, b);
		} while (a == b || h.find(make_pair(a, b)) != h.end());
		h.insert(make_pair(a, b));
		ret.push_back(make_pair(a, b));
	}
	return ret;
}

int main(int argc, char *argv[]) {
	if (argc != 5) {
		cout << "Usage: gen_test [filename] [n_terminals] [n_nets] [width/height]";
		return 0;
	}
	string filename = string(argv[1]);
	int N = atoi(argv[2]);
	int M = atoi(argv[3]);
	double cord = atof(argv[4]);

	ofstream fout(filename);
	auto point = random_points(N, cord);
	auto net = random_edges(N, M);
	fout << N << " " << M << " " << cord << " " << cord << endl;
	for (auto &x : point)
		fout << x.x << " " << x.y << endl;
	for (auto &x : net)
		fout << x.first << " " << x.second << endl;
	fout.close();

	return 0;
}