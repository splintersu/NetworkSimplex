#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;

#define mp make_pair 
#define fi first
#define se second

typedef int arc;

const int maxn = 500000;
const int maxm = 500000;

namespace solver
{
	const int special_cost = -1000000000;

	int n = -1 , m = 0;//node : [0 , n)
	//######### need initial value
	bool edge_from_sink_to_source_added = false;

	//arcs in residual net. [2 , 2m + 2)
	int nxt[maxm] , cap[maxm] , start_point[maxm] , to[maxm] , cost[maxm];
	int head_n[maxn] = {0};
	int head_t[maxn] = {0};
	bool tree_edge[maxm];
	int tot = 2;
	//####### need clear
	
	int total_cost = 0;
	int total_flow = 0;
	//####### need clear

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
		struct Edge{
			int st , ed , cost , capacity;
		}edges[maxm];
		
		int T_edges[maxn] , T_len = 0;
		
		int parent[maxn];
		int getf(int x){while(parent[x] != -1)x = parent[x];return x;}

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
			T_len = 0;
			for(int i = 0 ; i < n ; i++)
				remain[i] = 0;
			for(int i = 0 ; i < n ; i++)
				parent[i] = -1;
			for(int i = m ; i >= 1 ; i--)
			{
				int pu = getf(edges[i].st);
				int pv = getf(edges[i].ed);
				if(pu == pv) // L
				{
					x[i] = 0;
					addarc(edges[i].st , edges[i].ed , edges[i].capacity , edges[i].cost , false);
					addarc(edges[i].ed , edges[i].st , 0 , -edges[i].cost , false);
				}
				else
				{
					addtreearc(edges[i].st , edges[i].ed , i + 1000000000);
					addtreearc(edges[i].ed , edges[i].st , i);
					T_edges[++T_len] = i;
					parent[pu] = pv;
				}
			}
			//printf("n = %d, tree size %d\n" , n , T_len);

			dfs0(0 , -1);

			int e;
			total_cost = 0;
			for(int i = 1 ; i <= T_len ; i++)
			{
				e = T_edges[i];
				addarc(edges[e].st , edges[e].ed , edges[e].capacity - x[e], edges[e].cost , true);
				addarc(edges[e].ed , edges[e].st , x[e] , -edges[e].cost , true);
				total_cost += (int)x[e] * edges[e].cost;
			}
				
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

	arc arc_to_parent[maxn];int dep[maxn] , parent[maxn];

	pair<int,arc> bottlenect;

	void get_min_flow(int st , int ed)
	{
		while(dep[st] > dep[ed])
		{
			tree_arc_list[++tree_arc_list_top] = arc_to_parent[st];

			if(cap[arc_to_parent[st]] < bottlenect.fi)
			{
				bottlenect.fi = cap[arc_to_parent[st]];
				bottlenect.se = arc_to_parent[st];
			}
			st = parent[st];
		}
		while(dep[ed] > dep[st])
		{

			tree_arc_list[++tree_arc_list_top] = 1 ^ arc_to_parent[ed];

			if(cap[1 ^ arc_to_parent[ed]] < bottlenect.fi)
			{
				bottlenect.fi = cap[1 ^ arc_to_parent[ed]];
				bottlenect.se = 1 ^ arc_to_parent[ed];
			}
			ed = parent[ed];
		}
		while(st != ed)
		{
			tree_arc_list[++tree_arc_list_top] = arc_to_parent[st];
			tree_arc_list[++tree_arc_list_top] = 1 ^ arc_to_parent[ed];

			if(cap[arc_to_parent[st]] < bottlenect.fi)
			{
				bottlenect.fi = cap[arc_to_parent[st]];
				bottlenect.se = arc_to_parent[st];
			}
			st = parent[st];

			if(cap[1 ^ arc_to_parent[ed]] < bottlenect.fi)
			{
				bottlenect.fi = cap[1 ^ arc_to_parent[ed]];
				bottlenect.se = 1 ^ arc_to_parent[ed];
			}
			ed = parent[ed];
		}
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
		get_min_flow(to[non_tree_arc] , start_point[non_tree_arc]);
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
			dep[to[i]] = dep[o] + 1;
			parent[to[i]] = o;
			arc_to_parent[to[i]] = i ^ 1;
			pi[to[i]] = pi[o] - cost[i];
			dfs1(to[i] , o);
		}
	}

	arc check()
	{
		pi[0] = 0;
		dep[0] = 1;
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
		if(n == -1 || !edge_from_sink_to_source_added)
		{
			printf("n not intialized or the edge from sink to source is not added!\n");
			assert(false);
		}
		initial_solution::init();
		int count = 0;
		while(1)
		{
			arc non_tree_arc = check();
			if(non_tree_arc == -1)break;
			update(non_tree_arc);
			count++;
		}
		return mp(total_flow , total_cost);
	}

	void clear()
	{
		tot = 2;
		total_cost = 0;
		total_flow = 0;
		m = 0;
		for(int i = 0 ; i < n ; i++)
		{
			head_n[i] = 0;
			head_t[i] = 0;
		}
		solver::initial_solution::clear();
		n = -1;
		edge_from_sink_to_source_added = false;
	}

	void addedge(int from , int t , int cap , int wei)
	{
		if(wei == -1000000000)edge_from_sink_to_source_added = true;
		m++;
		initial_solution::edges[m].st = from;
		initial_solution::edges[m].ed = t;
		initial_solution::edges[m].capacity = cap;
		initial_solution::edges[m].cost = wei;
	}
}
