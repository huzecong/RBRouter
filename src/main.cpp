// For testing
#include <iostream>
#include <random>
#include <vector>
#include "RBRouter.h"

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

int main(int argc, char *argv[]) {
/*
	// ROAR
	vector<Point> vec = {{2.48, 2.56}, {6.6, 6.82},
						 {4.14, 10.14}, {5.84, 5.68},
						 {10.36, 5.26}, {4.18, 6.96}};
	RBRouter router(12.0, 12.0, vec);
	router.insert_net(0, 1);
	router.insert_net(2, 3);
	router.insert_net(4, 5);
*/
/*
	// normal (no border)
 	vector<Point> vec = {{1, 10}, {10, 10}, {20, 10},
						 {20, 20}, {10, 1}, {5, 5}, {15, 5}};
	RBRouter router(30.0, 30.0, vec);
	router.insert_net(0, 1);
	router.insert_net(1, 2);
	router.insert_net(2, 3);
	router.insert_net(1, 4);
	router.insert_net(5, 6);
*/
/*
	// Border
	vector<Point> vec = {{0, 10}, {10, 10}, {20, 10},
						 {20, 20}, {10, 0},
						 {5, 5}, {15, 5}};
	RBRouter router(20.0, 20.0, vec);
	router.insert_net(0, 1);
	router.insert_net(1, 2);
	router.insert_net(2, 3);
	router.insert_net(1, 4);
	router.insert_net(5, 6);
*/
/*
	// Rubber (thunder)
	vector<Point> vec = {{0, 10}, {20, 10},
						 {10, 10}, {30, 10},
						 {15, 0}, {15, 20}};
	RBRouter router(30.0, 20.0, vec);
	router.insert_net(0, 1);
	router.insert_net(2, 3);
	router.insert_net(4, 5);
*/
/*
	// Comprehensive
	vector<Point> vec = {{0, 100}, {5.15, 84.85},
						 {0, 27.06}, {0, 5.93},
						 {0, 15.83}, {10.62, 16.27},
						 {10, 100}, {56.15, 16.57},
						 {17.86, 100}, {49.64, 16.42},
						 {35.45, 100}, {45.5, 82.83},
						 {47.57, 91.21}, {30, 10},
						 {53.19, 100}, {39, 11.54}};
	RBRouter router(60.0, 100.0, vec);
	for (int i = 0; i < 8; ++i)
		router.insert_net(i * 2, i * 2 + 1);
*/
/*
	vector<Point> vec = {{11.36, 30}, {12.25, 0},
						 {0, 22.04}, {29.98, 21.6},
						 {21.26, 12.58}, {35.3, 12.14}};
	RBRouter router(40.0, 30.0, vec);
	for (int i = 0; i < 3; ++i)
		router.insert_net(i * 2, i * 2 + 1);
*/
/*
	vector<Point> vec = {{6.78, 8.44}, {19.05, 21.6}, {35.16, 21.6},
						 {29.1, 25.88}, {28.36, 13.76},
						 {5.15, 19.53}, {20, 10}};
	RBRouter router(40.0, 30.0, vec);
	router.insert_net(0, 1);
	router.insert_net(1, 2);
	router.insert_net(3, 4);
	router.insert_net(5, 6);
*/
/*
	vector<Point> vec = {{20, 60}, {44.97, 77.34}, {46.91, 21.74},
						 {75.68, 60.47}, {83.97, 79}, {72.08, 116.9}};
	RBRouter router(100.0, 140.0, vec);
	router.insert_net(1, 2);
	router.insert_net(3, 5);
	router.insert_net(0, 4);
	router.insert_net(1, 3);
*/

	const int N = 7;
	const double M = 50.0;
	vector<Point> vec = random_points(N, M);
	for (const Point &p : vec)
		debug(p.x << " " << p.y);
	RBRouter router(M, M, vec);
	for (int i = 1; i <= N * 1.5; ++i) {
		int a, b;
		do {
			a = rand() % N, b = rand() % N;
		} while (a == b);
		router.insert_net(a, b);
		debug(a << " " << b);
	}

	debug("Before router.solve()");
	double ans = router.solve();
	cout << ans << endl;

	router.plot("out.ps");

	return 0;
}