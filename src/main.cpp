// For testing
#include <iostream>
#include <vector>
#include "RBRouter.h"

using namespace std;

int main() {
/*
	vector<Point> vec = {{2.48, 2.56}, {6.6, 6.82},
						 {4.14, 10.14}, {5.84, 5.68},
						 {10.36, 5.26}, {5.18, 6.96}};
	RBRouter router(12.0, 12.0, vec);
	router.insert_net(0, 1);
	router.insert_net(2, 3);
	router.insert_net(4, 5);
*/
/*
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
	vector<Point> vec = {{0, 10}, {20, 10},
						 {10, 10}, {30, 10},
						 {15, 0}, {15, 20}};
	RBRouter router(30.0, 20.0, vec);
	router.insert_net(0, 1);
	router.insert_net(2, 3);
	router.insert_net(4, 5);
*/
/*
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
	vector<Point> vec = {{6.78, 8.44}, {19.05, 21.6}, {35.16, 21.6},
						 {29.1, 25.88}, {28.36, 13.76},
						 {5.15, 19.53}, {20, 10}};
	RBRouter router(40.0, 30.0, vec);
	router.insert_net(0, 1);
	router.insert_net(1, 2);
	router.insert_net(3, 4);
	router.insert_net(5, 6);

	debug("Before router.solve()");
	double ans = router.solve();
	cout << ans << endl;

	router.plot("out.ps");

	return 0;
}