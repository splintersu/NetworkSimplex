#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
#include "network_simplex.h"
using namespace std;

int u[100000] , v[100000] , c[100000] , w[100000];

int main()
{
	int n , m , source , sink;
	scanf("%d%d%d%d" , &n , &m , &source , &sink);
	solver::clear();
	solver::n = n;
	for(int i = 1 ; i <= m ; i++)
	{
		scanf("%d%d%d%d" , &u[i] , &v[i] , &c[i] , &w[i]);
		solver::addedge(u[i] , v[i] , c[i] , w[i]);
	}
	solver::addedge(sink , source , 1000000000 , -1000000000);
	auto res = solver::solve();
	printf("flow = %d, cost = %d\n" , res.first , res.second);
	/*
	solver::clear();
	solver::n = n;
	for(int i = m ; i >= 1 ; i--)
		solver::addedge(v[i] , u[i] , c[i] , w[i]);
	solver::addedge(source , sink , 1000000000 , -1000000000);
	res = solver::solve();
	printf("flow = %d, cost = %d\n" , res.fi , res.se);
	*/

	return 0;
}
