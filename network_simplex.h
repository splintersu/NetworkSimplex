#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;


namespace solver
{
	#define mp make_pair 
	#define fi first
	#define se second

	typedef int arc;

	const int maxn = 500000;
	const int maxm = 500000;


	const int special_cost = -1000000000;

	int n , m;//node : [0 , n)

	//arcs in residual net. [2 , 2m + 2)
	int nxt[maxm] , cap[maxm] , start_point[maxm] , to[maxm] , cost[maxm];
	int head_n[maxn] = {0};
	int head_t[maxn] = {0};
	bool tree_edge[maxm];int tot = 2;
	
	int total_cost;
	int total_flow;

	void addarc(int from , int t , int capacity , int weight ,  bool tree)
	{
		int now = tot++;
		tree_edge[now] = tree;
		to[now] = t;
		start_point[now] = from;
		cap[now] = capacity;
		cost[now] = weight;
		if(tree)
		{
			nxt[now] = head_t[from];
			head_t[from] = now;
		}
		else
		{
			nxt[now] = head_n[from];
			head_n[from] = now;
		}
	}

	namespace initial_solution
	{
		int b[maxn];
		//b stands for demanded flow, also the extra flow in
		struct Edge{
			int st , ed , cost , capacity;
		}edges[maxm];
		vector<int> L , U , T; // represents the initial feasible solution
		
		int parent[maxn];
		int getf(int x){while(parent[x] != -1)x = parent[x];return x;}
		void generate_random_spanning_tree()
		{
			for(int i = 0 ; i < n ; i++)
				parent[i] = -1;
			for(int i = m ; i >= 1 ; i--)
			{
				int pu = getf(edges[i].st);
				int pv = getf(edges[i].ed);
				if(pu == pv)L.push_back(i);
				else
				{
					T.push_back(i);
					parent[pu] = pv;
				}
			}
			//printf("n = %d, tree size %d\n" , n , int(T.size()));
		}

		int x[maxm];
		
		int head[maxn] = {0} , to[maxm] , nxt[maxm] , index[maxm] , tot = 0;
		
		void clear()
		{
			tot = 0;
			for(int i = 0 ; i < n ; i++)
				head[i] = 0;
		}

		void addtreearc(int from , int t , int ind)
		{
			int now = ++tot;
			nxt[now] = head[from];
			to[now] = t;
			index[now] = ind;
			head[from] = now;
		}

		int remain[maxn];//real flow in - real flow out

		void dfs0(int o , int fa)
		{
			for(int i = head[o] ; i ; i = nxt[i])
			{
				if(to[i] == fa)continue;
				dfs0(to[i] , o);
			}
			if(fa != -1)
			{
				int arcindex;
				for(int i = head[o] ; i ; i = nxt[i])
					if(to[i] == fa)
						arcindex = index[i];

				if(arcindex >= 1000000000) //forward arc, <o , fa>
				{
					//printf("forward<%d,%d>\n" , o , fa);
					arcindex -= 1000000000;
					x[arcindex] = remain[o];
					remain[o] -= x[arcindex];
					remain[fa] += x[arcindex];
				}
				else //backward arc, <fa , o>
				{
					//printf("backward<%d,%d>\n" , fa , o);
					x[arcindex] = -remain[o];
					remain[o] += x[arcindex];
					remain[fa] -= x[arcindex];
				}
			}
		}

		void init()
		{
			for(int i = 0 ; i < n ; i++)
				remain[i] = b[i];
			for(auto e : L)
				x[e] = 0;
			for(auto e : U)
			{
				x[e] = edges[e].capacity;
				remain[edges[e].ed] += edges[e].capacity;
				remain[edges[e].st] -= edges[e].capacity;
			}
			for(auto e : T)
			{
				addtreearc(edges[e].st , edges[e].ed , e + 1000000000);
				addtreearc(edges[e].ed , edges[e].st , e);
			}

			dfs0(0 , -1);
			/*
			for(int i = 1 ; i <= m ; i++)
				printf("x[%d] = %d\n" , i ,x[i]);
			*/
			for(auto e : T)
				if(x[e] < 0 || x[e] > edges[e].capacity)
					assert(false);

			for(auto e : L)
			{
				addarc(edges[e].st , edges[e].ed , edges[e].capacity , edges[e].cost , false);
				addarc(edges[e].ed , edges[e].st , 0 , -edges[e].cost , false);
			}
			for(auto e : U)
			{
				addarc(edges[e].st , edges[e].ed , 0 , edges[e].cost , false);
				addarc(edges[e].ed , edges[e].st , edges[e].capacity , -edges[e].cost , false);
			}
			
			for(auto e : T)
			{
				addarc(edges[e].st , edges[e].ed , edges[e].capacity - x[e], edges[e].cost , true);
				addarc(edges[e].ed , edges[e].st , x[e] , -edges[e].cost , true);
			}

			total_cost = 0;
			for(auto e : U)
				total_cost += (int)edges[e].capacity * edges[e].cost;
			for(auto e : T)
				total_cost += (int)x[e] * edges[e].cost;
			//printf("initial_cost = %lld\n" , total_cost);
		}
	}

	void remove_arc(int& head , int arc_index)
	{
		int first_arc = head;
		if(first_arc == arc_index)
		{
			head = nxt[arc_index];
			nxt[arc_index] = 0;
			return;
		}
		int pre_arc = first_arc;
		while(1)
		{
			first_arc = nxt[first_arc];
			if(first_arc == arc_index)
			{
				nxt[pre_arc] = nxt[arc_index];
				nxt[arc_index] = 0;
				return;
			}
			pre_arc = first_arc;
		}
	}
	void insert_arc(int& head , int arc_index)
	{
		nxt[arc_index] = head;
		head = arc_index;
	}

	int pi[maxn];

	int tree_arc_list[maxn] , tree_arc_list_top;

	pair<int,arc> bottlenect;

	bool get_min_flow(int o , int ed , int fa)
	{
		if(o == ed)return true;
		for(int i = head_t[o] ; i ; i = nxt[i])
		{
			if(to[i] == fa)continue;
			if(get_min_flow(to[i] , ed , o))
			{
				tree_arc_list[++tree_arc_list_top] = i;
				if(cap[i] < bottlenect.fi)
				{
					bottlenect.fi = cap[i];
					bottlenect.se = i;
				}
				return true;
			}
		}
		return false;
	}

	void augment(int f)
	{
		for(int i = 1 ; i <= tree_arc_list_top ; i++)
		{
			cap[tree_arc_list[i]] -= f;
			cap[tree_arc_list[i] ^ 1] += f;
			
			if(cost[tree_arc_list[i]] == special_cost)
					total_flow += f;
				else
					total_cost += (int)cost[tree_arc_list[i]] * f;
		}
		tree_arc_list_top = 0;
	}
	void update(arc non_tree_arc)
	{
		bottlenect = mp(2100000000 , -1);
		get_min_flow(to[non_tree_arc] , start_point[non_tree_arc] , -1);
		pair<int,arc> min_flow = bottlenect;
		//printf("bottlenect arc in tree <%d,%d>cap = %d\n" , min_flow.se.fi , to[min_flow.se.se] , min_flow.fi);
		if(cap[non_tree_arc] < min_flow.fi)
		{
			min_flow.fi = cap[non_tree_arc];
			min_flow.se = non_tree_arc;
		}
		//printf("bottlenect arc <%d,%d>cap = %d\n" , min_flow.se.fi , to[min_flow.se.se] , min_flow.fi);
		
		//augment(to[non_tree_arc] , start_point[non_tree_arc] , min_flow.fi , -1);
		augment(min_flow.fi);

		cap[non_tree_arc] -= min_flow.fi;
		cap[non_tree_arc ^ 1] += min_flow.fi;
		if(cost[non_tree_arc] == special_cost)
			total_flow += min_flow.fi;
		else
			total_cost += cost[non_tree_arc] * min_flow.fi;

		if(min_flow.se == non_tree_arc);
		else
		{
			remove_arc(head_n[start_point[non_tree_arc]] , non_tree_arc);
			remove_arc(head_n[to[non_tree_arc]] , non_tree_arc ^ 1);
			insert_arc(head_t[start_point[non_tree_arc]] , non_tree_arc);
			insert_arc(head_t[to[non_tree_arc]] , non_tree_arc ^ 1);

			arc remove_tree_arc = min_flow.se;
			remove_arc(head_t[start_point[remove_tree_arc]] , remove_tree_arc);
			remove_arc(head_t[to[remove_tree_arc]] , remove_tree_arc ^ 1);
			insert_arc(head_n[start_point[remove_tree_arc]] , remove_tree_arc);
			insert_arc(head_n[to[remove_tree_arc]] , remove_tree_arc ^ 1);

			tree_edge[non_tree_arc] = true;
			tree_edge[non_tree_arc ^ 1] = true;
			tree_edge[remove_tree_arc] = false;
			tree_edge[remove_tree_arc ^ 1] = false;
		}
	}

	void dfs1(int o , int fa)
	{
		for(int i = head_t[o] ; i ; i = nxt[i])
		{
			if(to[i] == fa)continue;
			pi[to[i]] = pi[o] - cost[i];
			dfs1(to[i] , o);
		}
	}

	arc check()
	{
		pi[0] = 0;
		dfs1(0 , -1);

		int min_c_pi = 0 , now_c_pi;
		arc min_arcindex;
		
		for(int i = 2 ; i < tot ; i++)
		{
			if(!tree_edge[i] && cap[i] > 0 && (now_c_pi = cost[i] - pi[start_point[i]] + pi[to[i]]) < min_c_pi)
			{
				min_c_pi = now_c_pi;
				min_arcindex = i;
			}
		}
		if(min_c_pi == 0)return -1;
		return min_arcindex;
	}

	pair<int , int> solve()
	{
		solver::initial_solution::init();
		int count = 0;
		while(1)
		{
			//printf("!\n");
			//calculate_pi();
			arc non_tree_arc = check();
			//printf("<%d,%d>\n" , start_point[non_tree_arc] , non_tree_arc);
			if(non_tree_arc == -1)break;
			//if(start_point[non_tree_arc] == -1)break;
			//printf("non_tree_arc is <%d,%d>, cap = %d\n" , start_point[non_tree_arc] , to[non_tree_arc] , cap[non_tree_arc]);
			update(non_tree_arc);
			count++;
		}
		return mp(total_flow , total_cost);
		//printf("total_cost : %lld\n" , total_cost);
	}

	void clear()
	{
		tot = 2;
		total_cost = 0;
		total_flow = 0;
		for(int i = 0 ; i < n ; i++)
		{
			head_n[i] = 0;
			head_t[i] = 0;
		}
		solver::initial_solution::clear();
	}

	void addedge(int from , int t , int cap , int wei)
	{
		solver::m++;
		solver::initial_solution::edges[solver::m].st = from;
		solver::initial_solution::edges[solver::m].ed = t;
		solver::initial_solution::edges[solver::m].capacity = cap;
		solver::initial_solution::edges[solver::m].cost = wei;
	}
}