#include <vector>
#include <map>
#include <optional>
#include <queue>
#include <tuple>
#include <functional>
#include <set>

using namespace std;

using Coord = pair<int, int>;
using Grid = vector<vector<int>>;

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

        for (const Coord& n_i : successor(n, t, grid)) {

            if (closed_nodes.count(n_i)) continue;
            double tentative_gscore = gscore[n] + 1.0;

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
