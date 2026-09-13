"""
PathPlanner.py

AERO60492 Coursework 2: Path Planning
A python module to implement the A* path planning algorithm with a Euclidean
distance heuristic function.

Author: 174347826+bt-nav@users.noreply.github.com
Date: March 2026

Using the pseudocode in [3], and the four steps outlined in [1], I have
implemented an A* algorithm in the do_a_star function. To improve the
algorithm given in [1], I have added the check tentative_g < g[n_i] to Step 4,
as suggested in [2]. This makes my implementation the "graph search" version of
A*, which is more efficient than the "tree search" version in [1], as it avoids
re-expanding nodes that have already been expanded with a cheaper path. To get
an optimal path from the graph search version, the heuristic must be admissible
and consistent [2], which is satisfied by Euclidean distance.

Python's min() function is used in Step 2 to find the node n with the smallest
f value. In the event of ties (i.e. f values are equally small), min() will
return the first node encountered. To ensure that the goal node t is favoured
in the event of ties, as specified in [1], I have used a lambda function as the
key for min(), which sorts by f[n] first, and then by n != t. For nodes that
are not t, n != t will evaluate to True, whereas t will evaluate to False. In
Python, False < True, so t will be favoured by min().

To optimise the algorithm further, one could focus on the data structure used
to store the open nodes. I have used a set, however nodes must be iterated over
until the node with the smallest f value is found. This results in a complexity
of O(n) for this step, where n is the number of open nodes. By using a
min-heap, the smallest f value can be retreived by popping the heap, reducing
complexity to a constant time O(1). I chose not to implement a min-heap as
importing libraries like heapq was not allowed for this coursework, and
implementation from scratch would have made the code overly complex.

Notation, from [1]:
n: a node
s: the start node
t: the goal node
T: the set of goal nodes, in this case T = {t}
Γ: the successor operator
g: the cost function,
   g(n) is the cost of the path from s to n
h: the heuristic function,
   h(n) is an estimate of the cost of an optimal path from n to t
f: the evaluation function,
   f(n) = g(n) + h(n), an estimate of the cost of a path from s to t through n.

References:
[1] P. Hart, N. Nilsson, and B. Raphael, 'A Formal Basis for the Heuristic
Determination of Minimum Cost Paths', IEEE Trans. Syst. Sci. Cyber., vol. 4,
no. 2, pp. 100-107, 1968, doi: 10.1109/TSSC.1968.300136.

[2] S. Russel and P. Norvig, 'A* search: Minimizing the total estimated
solution cost', in Artificial Intelligence: A Modern Approach, 3rd ed. Upper
Saddle River, NJ, 2009. Accessed: Mar. 20, 2026. [Online]. Available:
https://aima.cs.berkeley.edu/3rd-ed/

[3] 'A* search algorithm - Pseudocode', Wikipedia. Apr. 17, 2025. Accessed:
Mar. 20, 2026. [Online]. Available:
https://en.wikipedia.org/w/index.php?title=A*_search_algorithm&oldid=1340884788#Pseudocode
"""


def do_a_star(grid, start, end, display_message):
    """
    Executes the A* algorithm, following the four steps outlined in [1]

    1)  Mark s "open" and calculate f(s).
    2)  Select the open node n whose value of f is smallest. Resolve ties
        arbitrarily, but always in favor of any node n∈T.
    3)  If n∈T, mark n "closed" and terminate the algorithm.
    4)  Otherwise, mark n closed and apply the successor operator Γ to n.

        * Calculate distance to the successor n_i through n. If this tentative
        g value is less than the previously recorded g(n_i), or g(n_i) is not
        yet recorded, set this value as g(n_i) and proceed with Step 4.
        Otherwise, move on to the next successor as this one will not give a
        cheaper path so need not be expanded fully.

        Calculate f for each successor of n and mark as open each successor
        not already marked closed. Remark as open any closed node n_i, which is
        a successor of n and for which f(n_i) is smaller now than it was when
        n_i, was marked closed. Go to Step 2.

    * Step 4 has been modified according to [2], making this algorithm the
    graph search version of A*.

    Adapted from the pseudocode in [3].
    """
    s = start
    t = end
    path = []
    # As s is the start node, predecessor[s] = None.
    predecessor = {s: None}
    # The cost of the path from s to s is 0, so g[s] = 0.
    g = {s: 0}

    display_message(f"Planning a path from {s} to {t}", "INFO")

    # =============================== Step 1 ==================================
    # Mark s as open in the set of known nodes to be expanded.
    open_nodes = set([s])
    # Calculate f(s) = g(s) + h(s) = 0 + h(s) = h(s).
    f = {s: _heuristic(node=s, goal=t)}

    while len(open_nodes) > 0:
        # =========================== Step 2 ==================================
        # Select n with the smallest f(n). In a tie, t will be favoured.
        n = min(open_nodes, key=lambda n: (f[n], n != t))

        # =========================== Step 3 ==================================
        # If n is the goal node, reconstruct the path and terminate.
        if n == t:
            path = _reconstruct_path(predecessor, t)
            break

        # =========================== Step 4 ==================================
        # Otherwise, mark n as closed and apply the successor operator Γ.
        open_nodes.remove(n)
        for n_i in _successor(node=n, goal=t, grid=grid):
            # Each move to a successor has a cost of 1.
            tentative_g = g[n] + 1
            # Graph search check, see [2].
            if n_i not in g or tentative_g < g[n_i]:
                # Expand and record this path, it is the cheapest so far.
                predecessor[n_i] = n
                g[n_i] = tentative_g
                f[n_i] = g[n_i] + _heuristic(node=n_i, goal=t)
                # Open n_i. It is a successor of n and f(n_i) is smaller than
                # when n_i was marked closed. This meets the conditions in [1].
                if n_i not in open_nodes:
                    open_nodes.add(n_i)

    # ================================= End ===================================
    display_message("Planning completed.", "INFO")
    return path


# Helper functions


def _heuristic(node, goal):
    """
    The heuristic function. Calculates h(n), an estimate of the cost of an
    optimal path from n to a goal node t.

    In this case, the heuristic is the Euclidean distance from n to t.
    """
    n = node
    t = goal

    h = ((n[0] - t[0]) ** 2 + (n[1] - t[1]) ** 2) ** 0.5
    return h


def _reconstruct_path(predecessor, goal):
    """
    Reconstructs the path from the s to t using the predecessors.

    From [1] "We can then reconstruct a minimum cost path from s to t ...
    simply by chaining back from t to s through the pointers."
    """
    n = goal
    path = []

    # Chain backwards from t to s.
    while n is not None:
        path.append(n)
        n = predecessor[n]

    # Reverse to get the path from s to t.
    return path[::-1]


def _successor(node, goal, grid):
    """
    The successor operator Γ.

    Finds the successors (reachable neighbours) of a node in the grid.
    """
    successors = []
    col, row = node

    # Check up, down, left and right.
    if row > 0 and grid[col][row - 1] == 1:
        n_i = (col, row - 1)
        successors.append(n_i)
    if row < len(grid[0]) - 1 and grid[col][row + 1] == 1:
        n_i = (col, row + 1)
        successors.append(n_i)
    if col > 0 and grid[col - 1][row] == 1:
        n_i = (col - 1, row)
        successors.append(n_i)
    if col < len(grid) - 1 and grid[col + 1][row] == 1:
        n_i = (col + 1, row)
        successors.append(n_i)

    return successors
