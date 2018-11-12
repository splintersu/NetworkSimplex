/* This is an implementation of Network Simplex Algorithm
 * @Author: Jinyue Su(s.jinyue@gmail.com)
 * https://github.com/splintersu/NetworkSimplex
 * 
 * Network Simplex Algorithm is for Minimum-Cost-Flow problem
 * Problem Setting:
 * Denotes the network G = (N, A)
 * input: c_{i, j} for cost, u_{i, j} for capacity
 * the goal: min z(x) = sum_{(i, j) \in A} c_{i,j} \cdot x_{i,j}
 * satisfying:
 * 0 \leq x_{i, j} \leq u_{i, j} \forall {i, j} \in A
 * sum_{(i, j) \in A} x_{i, j} - sum_{(k, i) in \A} x_{k, i} = 0
 * 
 * Note that the problem setting is a little different from the one on wikipedia
 * This implementation follows
 * http://www.unc.edu/depts/stat-or/courses/provan/STOR724_web/lect14_simp.pdf
 * Note that, if there is no negative loop, it's possible to reduce Maximum-Flow-Minimum-Cost to MCFP:
 * just add an arc from sink to source, with infinity capacity and
 * -INF(some large enough value, take care not to cause overflow) weight
 * then call MCFP solver, decode the true flow and true cost from the returned cost
 */

#include <vector>

class NetworkSimplexSolver {
public:
	using cap_t = unsigned int; // type of capacity
	using cost_t = long long; // type of cost
	using pint_t = unsigned int; // type of positive integer

	struct Arc {
		pint_t st, ed;
		cap_t capacity;
		cost_t cost_per_unit;
	};

	/* @param n: the number of nodes
	 */

	/* @param n: number of nodes in the network
	 * @param arcs: Arcs in the network
	 * NOTE: The virtual reverse arc will be handled within this function,
	 * don't bother to add both arcs.
	 */
	static cost_t Solve(pint_t n, const std::vector<Arc>& arcs);
private:
	static bool SanityCheck(pint_t n, const std::vector<Arc>& arcs);
};
