/**
 * @file pathfinder.h
 * @author btnav (174347826+bt-nav@users.noreply.github.com)
 * @brief A C++ implementation of the A* pathfinding algorithm.
 *
 * Notation, from [1]:
 * n, a node.
 * s, the start node.
 * t, the target node.
 * T, the set of target nodes, in this case T = {t}.
 * Γ, the successor operator, which returns reachable neighbours.
 * g, the cost function.
 *    g(n) is the cost of the path from s to n.
 * h, the heuristic function.
 *    h(n) is the estimated cost of an optimal path from n to t.
 * f, the evaluation function.
 *    f(n) = g(n) + h(n), the estimated cost of a path from s to t, via n.
 *
 * References:
 * [1] P. Hart, N. Nilsson, and B. Raphael, 'A Formal Basis for the Heuristic
 * Determination of Minimum Cost Paths', IEEE Trans. Syst. Sci. Cyber., vol. 4,
 * no. 2, pp. 100-107, 1968, doi: 10.1109/TSSC.1968.300136.
 *
 * [2] S. Russel and P. Norvig, 'A* search: Minimizing the total estimated
 * solution cost', in Artificial Intelligence: A Modern Approach, 3rd ed. Upper
 * Saddle River, NJ, 2009. Accessed: Mar. 20, 2026. [Online]. Available:
 * https://aima.cs.berkeley.edu/3rd-ed/
 *
 * [3] 'A* search algorithm - Pseudocode', Wikipedia. Apr. 17, 2025. Accessed:
 * Mar. 20, 2026. [Online]. Available:
 * https://en.wikipedia.org/w/index.php?title=A*_search_algorithm&oldid=1340884788#Pseudocode
 *
 * @date September 2026
 *
 * @copyright Copyright (c) 2026, btnav (https://github.com/btnav).
 * All rights reserved. Released under the BSD 3-Clause License.
 *
 */
#ifndef PATHFINDER_H
#define PATHFINDER_H

#include <map>
#include <optional>
#include <tuple>
#include <vector>

using namespace std;

using Grid = vector<vector<int>>; // A 2D grid. 1 = empty, 0 = obstacle.
using Coord = pair<int, int>;  // A grid coordinate (col,row).
using Path = vector<Coord>; // A sequence of coordinates.

/**
 * @brief The heuristic function, h.
 *
 * Calculates h(n), the estimated cost of the optimal path from n to t.
 * Euclidean distance from n to t has been chosen as the heuristic.
 * 
 * @param node The current node, n.
 * @param target The target node, t.
 * @return double The Euclidean distance from n to t.
 */
double heuristic(const Coord& node, const Coord& target);

/**
 * @brief Finds a path from s to t.
 *
 * Reconstructs the minimum cost path from s to t by backtracking through
 * predecessors from t to s.
 *
 * @param prev Map of nodes to their predecessors.
 * @param target The target node, t.
 * @return Path Reconstructed path from s to t.
 */
Path reconstruct_path(
    const map<Coord, optional<Coord>>& prev,
    const Coord& target
);

/**
 * @brief The successor operator, Γ.

 * Finds the successors (reachable neighbours) of a node in the grid.

 * @param node The current node, n.
 * @param grid The grid to search.
 * @return Path Coordinates of n's successors.
 */
Path get_successors(
    const Coord& node,
    const Grid& grid
);

/**
 * @brief Executes the A* algorithm, following the steps outlined in [1].
 *
 * Adapted from the pseudocode in [3]:
 *
 * 1)  Mark s "open" and calculate f(s).
 * 2)  Select the open node n whose value of f is smallest. Resolve ties
 *     arbitrarily, but always in favor of any node n∈T.
 * 3)  If n∈T, mark n "closed" and terminate the algorithm.
 * 4)  Otherwise, mark n closed and apply the successor operator to n.
 *
 *     Calculate distance to the successor n_i through n. If this tentative
 *     g-score is less than the previously recorded g(n_i), or if g(n_i) is not
 *     yet recorded, set this value as g(n_i) and proceed with Step 4.
 *     Otherwise, move on to the next successor; this one will not give a
 *     cheaper path so need not be expanded fully. This is the graph search
 *     version of A*, outlined in [2].
 *
 *     Calculate f for each successor of n and mark as open each successor
 *     not already marked closed. Remark as open any closed node n_i, which is
 *     a successor of n and for which f(n_i) is smaller now than it was when
 *     n_i, was marked closed. Go to Step 2.
 *
 * @param grid The 2D grid to navigate.
 * @param start The starting node, s.
 * @param target The target node, t
 * @return Path The shortest path between s and t, if one exists.
 */
Path a_star(
    const Grid& grid,
    const Coord& start,
    const Coord& target
);

#endif // PATHFINDER_H