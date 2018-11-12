#include <bits/stdc++.h>
#include "network_simplex_solver.h"

using namespace std;

int main() {
	vector<NetworkSimplexSolver::Arc> arcs;
	arcs.push_back({0, 1, 1, -1});
	arcs.push_back({1, 2, 2, -1});
	arcs.push_back({2, 0, 3, -1});

	arcs.push_back({2, 3, 2, -5});
	arcs.push_back({3, 1, 2, 0});

	cout << NetworkSimplexSolver::Solve(4, arcs) << endl;
	return 0;
}