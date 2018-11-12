# Network Simplex
A C++ implementation of Network Simplex Algorithm

Author: splintersu@github

https://github.com/splintersu/NetworkSimplex

------

Network Simplex Algorithm is for Minimum-Cost-Flow problem

Problem Setting:

Denotes the network $G = (N, A)$

input: $c_{i, j}$ for cost, $u_{i, j}$ for capacity

the goal: min $z(x) = sum_{(i, j) \in A} c_{i,j} \cdot x_{i,j}$

satisfying:

$0 \leq x_{i, j} \leq u_{i, j} \forall {i, j} \in A$

$sum_{(i, j) \in A} x_{i, j} - sum_{(k, i) \in A} x_{k, i} = 0$

------

reference:

http://www.unc.edu/depts/stat-or/courses/provan/STOR724_web/lect14_simp.pdf