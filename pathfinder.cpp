/**
 * @file pathfinder.cpp
 * @author btnav (174347826+bt-nav@users.noreply.github.com)
 * @brief A C++ implementation of the A* pathfinding algorithm.
 * 
 * This A* pathfinding algorithm is based on the pseudocode in [3], and the
 * four steps outlined in [1]. To improve the algorithm given in [1], the check
 * tent_g < g[n_i] has been added to to Step 4, as suggested in [2]. This is the 
 * "graph search" version of A*, which is more efficient than the "tree search" 
 * version in [1], as it avoids re-expanding nodes that have already been 
 * expanded with a cheaper path. To get an optimal path from the graph search
 * version, the heuristic must be admissible and consistent [2], which is 
 * satisfied by Euclidean distance.
 *
 * A min-heap has been used to store the open nodes. When using a min-heap, the
 * smallest f-score can be accessed by topping the heap, with a constant time
 * complexity, O(1). Previously, a set was used. If using a set, nodes must be
 * iterated over until the node with the smallest f-score is found. This gives
 * a linear time complexity of O(n), which is less desirable.
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
#include "pathfinder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <tuple>
#include <vector>

using namespace std;

double heuristic(const Coord& node, const Coord& target) {
    int dx = (node.first - target.first);
    int dy = (node.second - target.second);

    // Calculate Euclidean distance.
    return static_cast<double>(sqrt(dx * dx + dy * dy));
}

Path reconstruct_path(
    const map<Coord, optional<Coord>>& prev,
    const Coord& target
) {
    Coord n = target;
    Path path;
    
    // Chain backwards from t to s.
    while (prev.at(n).has_value()) {
        path.push_back(n);
        n = *prev.at(n);
    }
    path.push_back(n);
    // Reverse to get the path from s to t.
    reverse(path.begin(), path.end());
    return path;
}

Path get_successors(
    const Coord& node,
    const Grid& grid
) {
    Path successors;
    int col = node.first;
    int row = node.second;

    // Check if up, down, left and right neigbours are accessible.
    if (row > 0 && grid[col][row - 1] == 1) {
        successors.emplace_back(col, row - 1);
    }
    if (row < (grid[0].size()) - 1 && grid[col][row + 1] == 1) {
        successors.emplace_back(col, row + 1);
    }
    if (col > 0 && grid[col - 1][row] == 1) {
        successors.emplace_back(col - 1, row);
    }
    if (col < (grid.size()) - 1 && grid[col + 1][row] == 1) {
        successors.emplace_back(col + 1, row);
    }
    return successors;
}

Path a_star(
    const Grid& grid,
    const Coord& start,
    const Coord& target
) {                      
    // Record the start time.
    auto start_time = chrono::steady_clock::now();
    
    // Initialise variables.
    Coord s = start; // the start node, s.
    Coord t = target; // the target node, t.
    Path path; // initial empty path.
    map<Coord, optional<Coord>> prev; // nodes to their predecessor.
    map<Coord, double> gscore; // nodes to their g-score.
    map<Coord, double> fscore; // nodes to their f-score.
    set<Coord> closed_nodes; // set of closed nodes.

    // Min-heap to store open nodes.
    using HeapEntry = tuple<double, int, Coord>;
    priority_queue<
        HeapEntry,
        vector<HeapEntry>,
        greater<HeapEntry>
    > open_nodes;
                  
    // The start node has no predecessor.
    prev[s] = nullopt;
    // The cost from s to s is 0, so g[s] = 0.
    gscore[s] = 0.0;

    // ============================== Step 1 ==================================
    
    // f(s) = g(s) + h(s) = 0 + h(s) = h(s).
    fscore[s] = heuristic(s, t);
    // Mark s as open.
    open_nodes.push({fscore[s], 0, s});

    while (!open_nodes.empty()) {
        // ========================== Step 2 ==================================
        
        // Select n with the smallest f-score. In a tie, t will be favoured.
        // Topping and popping the heap has constant time complexity O(1).
        auto [f, tie, n] = open_nodes.top();
        open_nodes.pop();

        // Skip closed or stale nodes.
        if (closed_nodes.count(n)) continue;
        if (fscore.count(n) && f > fscore[n]) continue;

        // ========================== Step 3 ==================================
        
        // If n is the goal node, reconstruct the path and terminate.
        if (n == t) {
            auto end_time = chrono::steady_clock::now();
            double elapsed_us = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            return reconstruct_path(prev, t);
        }

        // ========================== Step 4 ==================================
        
        // Otherwise, mark n as closed and apply the successor operator Γ.
        closed_nodes.insert(n);
        for (const Coord& n_i : get_successors(n, grid)) {

            // Skip closed nodes.
            if (closed_nodes.count(n_i)) continue;

            // Moving to a successor would have a cost of 1.
            double tent_gscore = gscore[n] + 1.0;
            // Graph search check, see [2].
            if (gscore.find(n_i) == gscore.end() || 
                tent_gscore < gscore[n_i]) {
                // Record this path, it is the cheapest so far.
                prev[n_i] = n;
                gscore[n_i] = tent_gscore;
                fscore[n_i] = tent_gscore + heuristic(n_i, t);
                // Open n_i, it meets the conditions in [1].
                open_nodes.push({fscore[n_i], (n_i == t) ? 0 : 1, n_i});
            }
        }

    }
    return path;
    // ================================ End ===================================
}

