#include <vector>
#include <map>
#include <optional>
#include <queue>
#include <tuple>
#include <functional>
#include <set>
#include <cmath>
#include <algorithm>

using namespace std;

using Coord = pair<int, int>;
using Grid = vector<vector<int>>;

double heuristic(const Coord& node, const Coord& goal) {
    int dx = (node.first - goal.first);
    int dy = (node.second - goal.second);
    return static_cast<double>(sqrt(dx * dx + dy * dy));
}

vector<Coord> reconstruct_path(
    const map<Coord, optional<Coord>>& prev,
    const Coord& goal
) {
    Coord n = goal;
    vector<Coord> path;
    
    while (prev.at(n).has_value()) {
        path.push_back(n);
        n = *prev.at(n);
    }
    path.push_back(n);
    
    reverse(path.begin(), path.end());
    return path;
}

vector<Coord> get_neighbours(
    const Coord& node,
    const Coord& goal,
    const Grid& grid
) {
    vector<Coord> neigbours;
    int col = node.first;
    int row = node.second;

    const int num_cols = static_cast<int>(grid.size());
    const int num_rows = static_cast<int>(grid[0].size());
    
    if (row > 0 && grid[col][row - 1] == 1) {
        neigbours.emplace_back(col, row - 1);
    }
    if (row < num_rows - 1 && grid[col][row + 1] == 1) {
        neigbours.emplace_back(col, row + 1);
    }
    if (col > 0 && grid[col - 1][row] == 1) {
        neigbours.emplace_back(col - 1, row);
    }
    if (col < num_cols - 1 && grid[col + 1][row] == 1) {
        neigbours.emplace_back(col + 1, row);
    }
    return neigbours;
}

vector<Coord> a_star(
    const Grid& grid,
    const Coord& start,
    const Coord& end
) {                      
    Coord s = start;
    Coord t = end;
    vector<Coord> path;
    map<Coord, optional<Coord>> prev;
    map<Coord, double> gscore;
    map<Coord, double> fscore;
    set<Coord> closed_nodes;
    
    using HeapEntry = tuple<double, int, Coord>;
    priority_queue<
        HeapEntry,
        vector<HeapEntry>,
        greater<HeapEntry>
    > open_nodes;
                  
    prev[s] = nullopt;
    gscore[s] = 0.0;
    fscore[s] = heuristic(s, t);
    open_nodes.push({fscore[s], (s == t) ? 0, s});

    while (!open_nodes.empty()) {
        auto [f, tie, n] = open_nodes.top();
        open_nodes.pop();

        if (closed_nodes.count(n)) continue;
        if (fscore.count(n) && f > fscore[n]) continue;

        if (n == t) {
            return reconstruct_path(prev, t);
        }

        closed_nodes.insert(n);

        for (const Coord& n_i : get_neighbours(n, t, grid)) {

            if (closed_nodes.count(n_i)) continue;
            double tent_gscore = gscore[n] + 1.0;

            if (gscore.find(n_i) == gscore.end() || 
                tent_gscore < gscore[n_i]) {

                prev[n_i] = n;
                gscore[n_i] = tent_gscore;
                fscore[n_i] = tent_gscore + heuristic(n_i, t);
                open_nodes.push({fscore[n_i], (n_i == t) ? 0 : 1, n_i});
            }
        }

    }
    return path;
}
