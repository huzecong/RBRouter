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

	vector<Point> vec = {{0, 10}, {20, 10},
						 {10, 10}, {30, 10},
						 {15, 0}, {15, 20}};
	RBRouter router(30.0, 20.0, vec);
	router.insert_net(0, 1);
	router.insert_net(2, 3);
	router.insert_net(4, 5);

	debug("Before router.solve()");
	double ans = router.solve();
	cout << ans << endl;

	router.plot("out.ps");

	return 0;
}