#include "network_simplex_solver.h"
#ifdef DEBUG
#include <cstdio>
#endif

using namespace std;
using cost_t = NetworkSimplexSolver::cost_t;
using cap_t = NetworkSimplexSolver::cap_t;
using pint_t = NetworkSimplexSolver::pint_t;

#ifdef DEBUG
#define LOG(args...) do{\
	printf("DEBUG: %s %d - <%s>: ", __FILE__, __LINE__, __FUNCTION__);\
	printf(args);printf("\n");\
}while(false);
#else
#define LOG(args...) do{}while(false);
#endif // DEBUG

struct Tree{
	vector<cost_t> pi;
	vector<pint_t> dep, parent, arc_to_parent;

	Tree(pint_t n) {
		pi = vector<cost_t>(n, 0);
		dep = parent = arc_to_parent = vector<pint_t>(n, 0);
	}

	void Clear() {
		dep[0] = 1;
		pi[0] = 0;
		parent[0] = -1;
		arc_to_parent[0] = -1;
	}

	template<class Func>
	void Traverse(pint_t st, pint_t ed, Func&& func){
		while(dep[st] > dep[ed]){
			func(arc_to_parent[st]);
			st = parent[st];
		}
		while(dep[st] < dep[ed]){
			func(arc_to_parent[ed] ^ 1);
			ed = parent[ed];
		}
		while(st != ed){
			func(arc_to_parent[st]);
			func(arc_to_parent[ed] ^ 1);
			st = parent[st];
			ed = parent[ed];
		}
	}
};

class Solution {
public:
	pint_t n;

	struct Arc {
		pint_t st, ed;
		cap_t capacity;
		cost_t cost_per_unit;

 		// These Arcs are stored in Linked-List.
 		// We cannot use pointer here since vector doesn't give persistent memory address.
		pint_t next, pre;

		// It's guaranteed by the caller that
		// (the index of an arc) xor (the index of the reverse arc) = 1

		Arc(NetworkSimplexSolver::Arc arc):
		st(arc.st), ed(arc.ed), capacity(arc.capacity), cost_per_unit(arc.cost_per_unit)
		{}

		Arc() = default;
	};
	vector<Arc> arcs;
	vector<pint_t> headn, headt;
	vector<bool> tree_edge;
	
	void AddArc(Arc arc, bool on_tree);
	void InsertArc(int arc_index, bool on_tree) {
		pint_t* head;
		if(on_tree)
			head = &headt[arcs[arc_index].st];
		else
			head = &headn[arcs[arc_index].st];

		arcs[arc_index].next = *head;
		arcs[arc_index].pre = 0;
		if(*head != 0)
			arcs[*head].pre = arc_index;
		*head = arc_index;
	}
	void RemoveArc(int arc_index, bool on_tree) {

		pint_t* head;
		if(on_tree)
			head = &headt[arcs[arc_index].st];
		else
			head = &headn[arcs[arc_index].st];

		int pre = arcs[arc_index].pre, next = arcs[arc_index].next;

		if(pre == 0)
			*head = next;
		else
			arcs[pre].next = next;
		
		if(next != 0)
			arcs[next].pre = pre;
	}

	cost_t cost;

	Solution(pint_t n): n(n), cost(0), tree(n) {
		headn = headt = vector<pint_t>(n, 0);
		arcs = vector<Arc>{Arc(), Arc()};
		tree_edge = vector<bool>{false, false};
		// to make sure the paired-arcs have indices: (2, 3), (4, 5), ...
	}

	static Solution GetZeroInitialSolution(pint_t n, const vector<NetworkSimplexSolver::Arc>& arcs);

	cost_t Solve();

private:

	Tree tree;
	// Tree information

	void Pivot(pint_t arc_to_augment) {
		pint_t min_capacity = arcs[arc_to_augment].capacity;
		pint_t argmin_capacity = arc_to_augment;

		auto get_min_capacity = [&](pint_t arc_index) {
			if(arcs[arc_index].capacity < min_capacity){
				min_capacity = arcs[arc_index].capacity;
				argmin_capacity = arc_index;
			}
		};
		tree.Traverse(arcs[arc_to_augment].ed, arcs[arc_to_augment].st, get_min_capacity);

		LOG("augment flow = %d, <%d,%d>", int(min_capacity),
			int(arcs[argmin_capacity].st), int(arcs[argmin_capacity].ed));
		cost_t cost_before_augment = cost;

		auto augment = [&](pint_t arc_index) {
			arcs[arc_index].capacity -= min_capacity;
			arcs[arc_index ^ 1].capacity += min_capacity;

			cost += arcs[arc_index].cost_per_unit * min_capacity;
		};
		tree.Traverse(arcs[arc_to_augment].ed, arcs[arc_to_augment].st, augment);
		arcs[arc_to_augment].capacity -= min_capacity;
		arcs[arc_to_augment ^ 1].capacity += min_capacity;
		cost += arcs[arc_to_augment].cost_per_unit * min_capacity;

		LOG("cost += %d", int(cost - cost_before_augment));

		if(argmin_capacity != arc_to_augment){
			RemoveArc(arc_to_augment, false);
			RemoveArc(arc_to_augment ^ 1, false);
			InsertArc(arc_to_augment, true);
			InsertArc(arc_to_augment ^ 1, true);

			pint_t remove_tree_arc = argmin_capacity;
			RemoveArc(remove_tree_arc, true);
			RemoveArc(remove_tree_arc ^ 1, true);
			InsertArc(remove_tree_arc, false);
			InsertArc(remove_tree_arc ^ 1, false);

			tree_edge[arc_to_augment] = true;
			tree_edge[arc_to_augment ^ 1] = true;
			tree_edge[remove_tree_arc] = false;
			tree_edge[remove_tree_arc ^ 1] = false;
		}
	}

	void DfsBuildTree(pint_t o, pint_t pa, Tree& tree) {
		for(pint_t i = headt[o] ; i != 0 ; i = arcs[i].next) {
			if(arcs[i].ed == pa)continue;

			int ed = arcs[i].ed;
			tree.dep[ed] = tree.dep[o] + 1;
			tree.parent[ed] = o;
			tree.arc_to_parent[ed] = i ^ 1;
			tree.pi[ed] = tree.pi[o] - arcs[i].cost_per_unit;

			DfsBuildTree(ed, o, tree);
		}
	}

	void BuildTree(){
		tree.Clear();
		DfsBuildTree(1, -1, tree);
	}

	pint_t FindArcToAugment(){

		cost_t min_c_pi = 0, now_c_pi;
		pint_t min_arcindex;

		for(pint_t i = 2 ; i < arcs.size() ; i++){
			if(!tree_edge[i] && arcs[i].capacity > 0){
				now_c_pi = arcs[i].cost_per_unit - tree.pi[arcs[i].st] + tree.pi[arcs[i].ed];
				if(now_c_pi < min_c_pi){
					min_c_pi = now_c_pi;
					min_arcindex = i;
				}
			}
		}

		if(min_c_pi == 0)
			// Cannot find
			return 0;
		return min_arcindex;
	}
};

void Solution::AddArc(Arc arc, bool on_tree){
	arcs.push_back(arc);
	tree_edge.push_back(on_tree);
	pint_t added_arc_index = arcs.size() - 1;

	InsertArc(added_arc_index, on_tree);
}

cost_t Solution::Solve() {
	while(true){
		BuildTree();
		pint_t arc_to_augment = FindArcToAugment();

		LOG("arc_to_augment = %d <%d,%d>", 
			int(arc_to_augment), int(arcs[arc_to_augment].st), int(arcs[arc_to_augment].ed));

		if(arc_to_augment == 0)break;
		Pivot(arc_to_augment);
	}

	return cost;
}

namespace GetInitialZeroSolutionDetails {

	class UnionFindSet{
	public:
		UnionFindSet(pint_t n): n(n) {
			p = vector<pint_t>(n, -1);
		}

		pint_t Find(pint_t a){
			if(p[a] == -1)
				return a;
			else
				return p[a] = Find(p[a]);
		}

		bool Union(pint_t a, pint_t b){
			pint_t pa = Find(a), pb = Find(b);
			if(pa == pb)return false;
			p[pa] = pb;
			return true;
		}

	private:
		pint_t n;
		vector<pint_t> p;
	};

	static Solution GetInitialSolution(pint_t n, const vector<NetworkSimplexSolver::Arc>& arcs) {
		Solution sol(n);

		LOG("call GetInitialSolution");

		UnionFindSet ufs(n);
		for(pint_t i = 0 ; i < arcs.size() ; i++){
			NetworkSimplexSolver::Arc reverse_arc = arcs[i];
			swap(reverse_arc.st, reverse_arc.ed);
			reverse_arc.capacity = 0;
			reverse_arc.cost_per_unit = -arcs[i].cost_per_unit;

			bool on_tree = ufs.Union(arcs[i].st, arcs[i].ed);
			sol.AddArc(arcs[i], on_tree);
			sol.AddArc(reverse_arc, on_tree);

			LOG("<%d,%d> capacity=%d, cost=%d, on the tree=%d",
				int(arcs[i].st), int(arcs[i].ed), int(arcs[i].capacity), int(arcs[i].cost_per_unit),
				int(on_tree));
		}

		return sol;
	}

}

// static
Solution Solution::GetZeroInitialSolution(pint_t n, const vector<NetworkSimplexSolver::Arc>& arcs) {
	return GetInitialZeroSolutionDetails::GetInitialSolution(n, arcs);
}

// static
cost_t NetworkSimplexSolver::Solve(pint_t n, const vector<Arc>& arcs) {
	Solution sol = Solution::GetZeroInitialSolution(n, arcs);
	LOG("finish build solution");
	return sol.Solve();
}
