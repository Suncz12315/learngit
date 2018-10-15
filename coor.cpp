#include"stdafx.h"
#include"coor.h"
vector<int>random_coor() {
	vector<int>coor(2,0);
	int a = 0, b = 5000;
	coor[0] = (rand() % (b - a)) + a;
	coor[1] = (rand() % (b - a)) + a;
	return coor;
}