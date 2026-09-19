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
    return static_cast<double>(sqrt(dx * dx + dy * dy));
}

Path reconstruct_path(
    const map<Coord, optional<Coord>>& prev,
    const Coord& target
) {
    Coord n = target;
    Path path;
    
    while (prev.at(n).has_value()) {
        path.push_back(n);
        n = *prev.at(n);
    }
    path.push_back(n);
    
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

    const int num_cols = static_cast<int>(grid.size());
    const int num_rows = static_cast<int>(grid[0].size());
    
    if (row > 0 && grid[col][row - 1] == 1) {
        successors.emplace_back(col, row - 1);
    }
    if (row < num_rows - 1 && grid[col][row + 1] == 1) {
        successors.emplace_back(col, row + 1);
    }
    if (col > 0 && grid[col - 1][row] == 1) {
        successors.emplace_back(col - 1, row);
    }
    if (col < num_cols - 1 && grid[col + 1][row] == 1) {
        successors.emplace_back(col + 1, row);
    }
    return successors;
}

Path a_star(
    const Grid& grid,
    const Coord& start,
    const Coord& target
) {                      
    auto start_time = chrono::steady_clock::now();
    
    Coord s = start;
    Coord t = target;
    Path path;
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
    open_nodes.push({fscore[s], 0, s});

    while (!open_nodes.empty()) {
        auto [f, tie, n] = open_nodes.top();
        open_nodes.pop();

        if (closed_nodes.count(n)) continue;
        if (fscore.count(n) && f > fscore[n]) continue;

        if (n == t) {
            auto end_time = chrono::steady_clock::now();
            double elapsed_us = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            return reconstruct_path(prev, t);
        }

        closed_nodes.insert(n);

        for (const Coord& n_i : get_successors(n, grid)) {

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

