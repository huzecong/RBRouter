//
// Created by Kanari on 15/7/3.
//

#include <vector>
#include <map>
#include "RBBase.h"

using namespace std;

RBNet RBNet::subnet(vector<ID> vec) const {
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
	vector<Point> sub_point;
	map<ID, ID> has_ID;
	for (unsigned int i = 0; i < vec.size(); ++i) {
		sub_point.push_back(this->point[vec[i]]);
		has_ID[vec[i]] = i;
	}
	RBNet result(this->width(), this->height(), sub_point);
	for (unsigned int i = 0; i < vec.size(); ++i)
		for (ID j : this->link[vec[i]]) {
			map<ID, ID>::iterator it = has_ID.find(j);
			if (it != has_ID.end() && it->second > i)
				result.add_net(i, it->second);
		}
	return result;
}