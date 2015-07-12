// For testing
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>
#include "RBRouter.h"

using namespace std;

int main(int argc, char *argv[]) {

	if (argc != 2 && argc != 3) {
		cout << "Usage: RBRouter [input_file] [output_postscript_file]" << endl;
		cout << "By default, [output_postscript_file] is \"out.ps\"" << endl;
		return 0;
	}
	string input_file = string(argv[1]);
	string output_file;
	if (argc == 3) output_file = string(argv[2]);
	else output_file = "out.ps";

	ifstream fin(input_file);
	int N, M;
	double width, height;
	fin >> N >> M >> width >> height;
	vector<Point> vec;
	for (int i = 0; i < N; ++i) {
		double x, y;
		fin >> x >> y;
		vec.push_back(Point(x, y));
	}
	RBRouter router(width, height, vec);
	for (int i = 0; i < M; ++i) {
		int x, y;
		fin >> x >> y;
		router.insert_net(x, y);
	}
	fin.close();

/*
	const int N = 4;
	const double M = 40.0;
	vector<Point> vec = random_points(N, M);
	for (const Point &p : vec)
		debug(p.x << " " << p.y);
	RBRouter router(M, M, vec);
	auto edges = random_edges(N, N * 1.5);
	for (auto x : edges) {
		router.insert_net(x.first, x.second);
		debug(x.first << " " << x.second);
	}
*/
	debug("Before router.solve()");
	double ans = router.solve();
	cout << ans << endl;

	if (ans != INFI) router.plot(output_file);

	return 0;
}