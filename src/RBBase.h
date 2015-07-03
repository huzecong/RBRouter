//
// Created by Kanari on 15/7/3.
//

#ifndef RBROUTER_RBMATH_H
#define RBROUTER_RBMATH_H


#include <cstdio>

struct Point {
	double x, y;
	Point() {}
	Point(double _x, double _y) : x(_x), y(_y) {}
};

struct Segment {
	Point s, e;
};

void ensure(const bool cond, const char *msg) {
	if (!cond) perror(msg);
}

#endif //RBROUTER_RBMATH_H
