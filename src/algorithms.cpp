#include <algorithm>
#include <algorithms.hpp>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <stdio.h>

namespace algorithms {

void get_isolated_vertices(Graph &graph, std::vector<int> &isolated_vertices) {
    isolated_vertices.clear();

    for (int v = 0; v < graph.size(); v++) {
        if (graph.valid_vertex(v)) {
            auto deg = graph.count_neighbors(v);

            if (deg == 0) {
                isolated_vertices.push_back(v);
            }
        }
    }
}

int clique_cover(const Graph &graph, std::vector<std::vector<int>> &cliques) {
    cliques.clear();

    int graph_size = graph.size();
    std::vector<int> clique_id(graph_size);

    for (int i = 0; i < graph_size; i++) {
        if (!graph.valid_vertex(i))
            continue;

        bool found_clique = false;
        auto adjacency = graph.get_neighbors(i);

        for (auto j : adjacency) {
            if (j > i)
                break;

            auto &clique = cliques[clique_id[j]];
            bool is_part_of_clique = true;

            for (const auto &v : clique) {
                if (!graph.has_neighbor(i, v)) {
                    is_part_of_clique = false;
                    break;
                }
            }

            if (is_part_of_clique) {
                clique.push_back(i);
                clique_id[i] = &clique - &(*cliques.begin());
                found_clique = true;
                break;
            }
        }

        if (!found_clique) {
            clique_id[i] = cliques.size();
            cliques.push_back({i});
        }
    }

    return graph.num_of_vertices() - cliques.size();
}

void find_satellites(const Graph &graph, int v, std::vector<int> &satellites) {
    auto adjacency_v = graph.get_neighbors(v);

    satellites.clear();

    for (const auto &w : adjacency_v) {
        if (!graph.valid_vertex(w))
            continue;

        auto adjacency_w = graph.get_neighbors(w);

        auto it = adjacency_v.begin();
        int u = INT_MAX;

        bool is_satellite = true;

        for (auto satellite : adjacency_w) {
            if (satellite == v)
                continue;

            while (it != adjacency_v.end() && *it < satellite)
                it++;

            if (it != adjacency_v.end() && *it == satellite)
                continue;

            if (u == INT_MAX) {
                u = satellite;
                continue;
            }

            is_satellite = false;
        }

        if (!is_satellite || u == INT_MAX)
            continue;

        satellites.push_back(u);
    }

    // Remove duplicates.
    std::sort(satellites.begin(), satellites.end());
    satellites.erase(unique(satellites.begin(), satellites.end()), satellites.end());
}

bool check_clique(const Graph &graph, std::vector<int> &clique) {
    for (const auto &x : clique)
        assert(graph.valid_vertex(x));

    for (auto x : clique) {
        for (auto y : clique) {
            if (x == y)
                continue;
            if (!graph.has_neighbor(x, y))
                return false;
        }
    }

    return true;
}

void find_mirrors(const Graph &graph, int v, std::vector<int> &mirrors) {
    mirrors.clear();

    std::vector<char> marked_v;
    std::vector<char> marked_u;
    marked_v.resize(graph.size(), 0);
    marked_u.resize(graph.size(), 0);

    auto adjacency_v = graph.get_neighbors(v);

    marked_v[v] = 1;
    for (auto w : adjacency_v) {
        marked_v[w] = 1;
    }

    auto check_mirror = [&]() -> bool {
        int n = 0;
        for (auto w : adjacency_v) {
            if (marked_u[w])
                continue;
            ++n;
        }
        for (auto w : adjacency_v) {
            if (marked_u[w])
                continue;
            int nx = 0;
            for (auto x : graph.get_neighbors(w)) {
                if (!marked_v[x] || marked_u[x])
                    continue;
                ++nx;
            }
            if (nx != n)
                return false;
        }
        return true;
    };

    for (auto w : adjacency_v) {
        for (auto u : graph.get_neighbors(w)) {
            if (u == v)
                continue;
            if (graph.has_neighbor(v, u))
                continue;

            auto adjacency_u = graph.get_neighbors(u);
            for (auto x : adjacency_u)
                marked_u[x] = 1;

            if (check_mirror())
                mirrors.push_back(u);

            for (auto x : adjacency_u)
                marked_u[x] = 0;
        }
    }

    // Remove duplicates.
    std::sort(mirrors.begin(), mirrors.end());
    mirrors.erase(unique(mirrors.begin(), mirrors.end()), mirrors.end());
}

void find_mirrors_and_satellites(const Graph &graph, int v, std::vector<int> &mirrors,
                                 std::vector<int> &satellites) {
    mirrors.clear();

    auto adjacency_v = graph.get_neighbors(v);

    for (const auto &w : adjacency_v) {
        auto adjacency_w = graph.get_neighbors(w);

        auto it = adjacency_v.begin();
        int satellite_id = INT_MAX;
        bool is_satellite = true;

        for (const auto &u : adjacency_w) {
            if (u == v)
                continue;

            if (!graph.has_neighbor(v, u)) {
                auto adjacency_u = graph.get_neighbors(u);

                std::vector<int> clique;
                std::set_difference(adjacency_v.begin(), adjacency_v.end(), adjacency_u.begin(),
                                    adjacency_u.end(), std::back_inserter(clique));

                if (check_clique(graph, clique))
                    mirrors.push_back(u);
            }

            // Satellites.
            while (it != adjacency_v.end() && *it < u)
                it++;

            if (it != adjacency_v.end() && *it == u)
                continue;

            if (satellite_id == INT_MAX) {
                satellite_id = u;
                continue;
            }

            is_satellite = false;
        }

        if (!is_satellite || satellite_id == INT_MAX)
            continue;

        satellites.push_back(satellite_id);
    }

    // Remove duplicates.
    std::sort(mirrors.begin(), mirrors.end());
    mirrors.erase(unique(mirrors.begin(), mirrors.end()), mirrors.end());
}

std::vector<int> upper_bound_2_guarantee(const Graph &graph) {
    int graph_size = graph.size();
    std::vector<int> in_cover(graph_size);
    std::vector<int> vertex_cover;

    for (int i = 0; i < graph_size; i++) {
        if (graph.valid_vertex(i) && !in_cover[i]) {
            auto neighbors = graph.get_neighbors(i);

            for (auto x : neighbors) {
                if (!in_cover[x]) {
                    in_cover[x] = 1;
                    in_cover[i] = 1;
                    vertex_cover.push_back(x);
                    vertex_cover.push_back(i);
                    break;
                }
            }
        }
    }

    return vertex_cover;
}

std::vector<int> upper_bound_greedy(const Graph &graph) {
    std::vector<int> vertices(graph.size());

    iota(vertices.begin(), vertices.end(), 0);
    vertices.erase(copy_if(vertices.begin(), vertices.end(), vertices.begin(),
                           [&](int x) { return graph.valid_vertex(x); }),
                   vertices.end());

    std::sort(vertices.begin(), vertices.end(), [&](int x, int y) {
        auto nx = graph.count_neighbors(x);
        auto ny = graph.count_neighbors(y);
        return nx > ny;
    });

    std::vector<int> in_cover(graph.size());
    std::vector<int> vertex_cover;

    for (const auto &x : vertices) {
        auto neighbors = graph.get_neighbors(x);

        for (auto y : neighbors) {
            if (!in_cover[y]) {
                in_cover[x] = 1;
                vertex_cover.push_back(x);
                break;
            }
        }
    }

    return vertex_cover;
}

std::vector<vertex_type> iterated_local_search(const Graph &g) {
    std::vector<int> state;
    std::vector<int> tightness; // Number of neighbors not in the VC.
    state.resize(g.size());
    tightness.resize(g.size());
    int n = 0;

    auto clear = [&] {
        n = 0;
        for (vertex_type u : g.get_valid_vertices()) {
            state[u] = 0;
            tightness[u] = g.count_neighbors(u);
        }
    };

    auto insert = [&](vertex_type x) {
        assert(!state[x]);
        state[x] = 1;
        ++n;

        for (vertex_type v : g.get_neighbors(x))
            --tightness[v];
    };

    auto remove = [&](vertex_type x) {
        assert(state[x]);
        state[x] = 0;
        --n;

        for (vertex_type v : g.get_neighbors(x))
            ++tightness[v];
    };

    // We start from a greedy solution. This is arbitrary.
    for (vertex_type u : upper_bound_greedy(g))
        insert(u);

    auto improve = [&](vertex_type x) -> bool {
        assert(!state[x]);

        for (vertex_type v : g.get_neighbors(x)) {
            if (!state[v])
                continue;
            if (tightness[v] != 1)
                continue;
            for (vertex_type w : g.get_neighbors(x)) {
                if (v == w)
                    continue;
                if (!state[w])
                    continue;
                if (tightness[w] != 1)
                    continue;
                if (g.has_neighbor(v, w))
                    continue;
                insert(x);
                remove(v);
                remove(w);
                return true;
            }
        }
        return false;
    };

    std::cerr << "before LS: " << n << std::endl;

    auto local_search = [&] {
        bool success;
        do {
            success = false;
            for (vertex_type x : g.get_valid_vertices()) {
                if (state[x])
                    continue;
                if (improve(x))
                    success = true;
            }
        } while (success);
    };

    local_search();
    std::vector<int> best_vc;
    for (vertex_type u : g.get_valid_vertices()) {
        if (state[u])
            best_vc.push_back(u);
    }

    std::cerr << "initial LS yields: " << best_vc.size() << std::endl;

    std::mt19937 prng(42);
    std::geometric_distribution<> pertub_distrib;
    std::uniform_int_distribution<> vertex_distrib(0, g.size() - 1);
    int i = 0;
    int i_success = 0;
    while (i < i_success + 50) {
        // Revert back to the best known VC.
        clear();
        for (vertex_type u : best_vc)
            insert(u);

        int m = pertub_distrib(prng) + 1;
        for (int j = 0; j < m; ++j) {
            vertex_type x;
            do {
                x = vertex_distrib(prng);
            } while (!g.valid_vertex(x) || !state[x]);

            remove(x);
            for (vertex_type v : g.get_neighbors(x)) {
                if (!state[v])
                    insert(v);
            }
        }

        // Remove non-tight vertices.
        for (vertex_type u : g.get_valid_vertices()) {
            if (!state[u])
                continue;
            if (!tightness[u])
                remove(u);
        }

        local_search();

        // If we find an improvement: capture the newly found VC.
        if (n < static_cast<int>(best_vc.size())) {
            best_vc.clear();
            for (vertex_type u = 0; u < g.size(); u++) {
                if (state[u])
                    best_vc.push_back(u);
            }
            i_success = i;
        }

        ++i;
    }

    std::cerr << "iterated LS yields: " << best_vc.size() << " after " << i << " iterations"
              << std::endl;

    return best_vc;
}

void get_connected_components(const Graph &graph,
                              std::vector<std::vector<int>> &connected_components) {
    int graph_size = graph.size();
    std::vector<int> visited(graph_size);
    std::queue<int> q;
    connected_components.clear();

    for (int i = 0; i < graph_size; i++) {
        if (!graph.valid_vertex(i) || visited[i])
            continue;

        q.push(i);
        visited[i] = true;
        connected_components.push_back({i});

        while (!q.empty()) {
            auto x = q.front();
            q.pop();
            auto neighbors = graph.get_neighbors(x);

            for (auto y : neighbors) {
                if (visited[y])
                    continue;

                visited[y] = true;
                connected_components.back().push_back(y);
                q.push(y);
            }
        }

        std::sort(connected_components.back().begin(), connected_components.back().end());
    }
}

bool dfs(int x, bool side, const Graph &G, std::vector<int> &matching,
         std::vector<int> &revmatching, const std::vector<int> &dist,
         const std::vector<int> &revdist, std::vector<int> &vis, std::vector<int> &revvis) {
    if (side == 0 && dist[x] == 0) {
        vis[x] = 1;
        return true;
    }
    if (side) {
        revvis[x] = 1;
    } else {
        vis[x] = 1;
    }

    if (side) {
        auto neighbors = G.get_neighbors(x);
        for (const auto &pre : neighbors) {
            if (vis[pre])
                continue;
            if (dist[pre] == revdist[x] - 1) {
                if (dfs(pre, 0, G, matching, revmatching, dist, revdist, vis, revvis)) {
                    matching[pre] = x;
                    revmatching[x] = pre;
                    return true;
                }
            }
        }
    } else {
        int pre = matching[x];
        if (revvis[pre])
            return false;
        if (revdist[pre] == dist[x] - 1) {
            if (dfs(pre, 1, G, matching, revmatching, dist, revdist, vis, revvis)) {
                return true;
            }
        }
    }
    return false;
}

void hopcroftKarp(const Graph &G, std::vector<int> &matching, std::vector<int> &revmatching) {
    int n = G.size();
    matching.assign(n, INT_MAX);
    revmatching.assign(n, INT_MAX);
    std::vector<int> dist(n, INT_MAX), revdist(n, INT_MAX);
    int sap = INT_MAX; // length of shortest augmenting path
    std::queue<int> q;
    std::vector<int> vis(n), revvis(n);

    do {
        fill(dist.begin(), dist.end(), INT_MAX);
        fill(revdist.begin(), revdist.end(), INT_MAX);

        sap = INT_MAX;
        while (!q.empty())
            q.pop();

        for (int i = 0; i < n; i++) {
            if (!G.valid_vertex(i))
                continue;
            if (G.valid_vertex(i) && matching[i] == INT_MAX) {
                q.push(i);
                dist[i] = 0;
            }
        }

        while (!q.empty()) {
            int t = q.front();
            q.pop();
            if (dist[t] + 1 > sap) {
                break;
            }
            auto neighbors = G.get_neighbors(t);
            for (const auto &nx : neighbors) {
                if (revdist[nx] != INT_MAX)
                    continue;
                if (nx == matching[t])
                    continue;

                revdist[nx] = dist[t] + 1;

                if (revmatching[nx] != INT_MAX) {
                    int t2 = revmatching[nx];
                    assert(t2 != t);

                    dist[t2] = dist[t] + 2;
                    q.push(t2);
                } else {
                    sap = dist[t] + 1;
                }
            }
        }
        if (sap != INT_MAX) {
            assert(sap & 1);
            fill(vis.begin(), vis.end(), 0);
            fill(revvis.begin(), revvis.end(), 0);
            bool aug = false;
            for (int i = 0; i < n; i++) {
                if (revdist[i] != INT_MAX && revmatching[i] == INT_MAX) {
                    assert(revdist[i] == sap);
                    if (dfs(i, 1, G, matching, revmatching, dist, revdist, vis, revvis)) {
                        aug = true;
                    }
                }
            }
            assert(aug);
        }
    } while (sap != INT_MAX);
}

void buildCover(std::vector<int> &cover, std::vector<int> &revcover, const Graph &G,
                const std::vector<int> &matching, const std::vector<int> &revmatching) {
    int n = G.size();
    cover.assign(n, 0);
    revcover.assign(n, 0);
    std::vector<int> Z(n), revZ(n);
    std::queue<int> q;
    for (int i = 0; i < n; i++) {
        if (G.valid_vertex(i) && matching[i] == INT_MAX) {
            q.push(i);
            Z[i] = 1;
        }
    }
    while (!q.empty()) {
        int t = q.front();
        q.pop();
        auto neighbors = G.get_neighbors(t);
        for (auto nx : neighbors) {
            if (nx == matching[t]) {
                continue;
            }
            int t2 = revmatching[nx];
            if (revZ[nx])
                continue;
            revZ[nx] = 1;
            Z[t2] = 1;
            q.push(t2);
        }
    }
    for (int i = 0; i < n; i++) {
        if (!G.valid_vertex(i))
            continue;
        cover[i] = 1 - Z[i];
        revcover[i] = revZ[i];
    }
}

int cycleCoverILPRelaxation(const Graph &G) {
    std::vector<int> matching;
    std::vector<int> revmatching;
    std::vector<int> cover;
    std::vector<int> revcover;

    hopcroftKarp(G, matching, revmatching);
    buildCover(cover, revcover, G, matching, revmatching);

    int n = G.size();
    int cycleBound = 0;
    std::vector<int> visited(G.size());
    for (int i = 0; i < n; i++) {
        if (G.valid_vertex(i) && !visited[i] && cover[i] + revcover[i] == 1) {
            assert(matching[i] != INT_MAX && revmatching[i] != INT_MAX);
            int length = 0;
            int x = i;
            do {
                visited[x] = 1;
                x = matching[x];
                length++;
            } while (x != i);
            cycleBound += length / 2;
        }
    }
    return cycleBound;
}

int cycleCoverILPRelaxation(const Graph &G, std::vector<int> &add, std::vector<int> &remove) {
    add.clear();
    remove.clear();

    std::vector<int> matching;
    std::vector<int> revmatching;
    std::vector<int> cover;
    std::vector<int> revcover;

    hopcroftKarp(G, matching, revmatching);
    buildCover(cover, revcover, G, matching, revmatching);

    int n = G.size();
    for (int i = 0; i < n; i++) {
        if (cover[i] && revcover[i])
            add.push_back(i);
        else if (G.valid_vertex(i) && !cover[i] && !revcover[i])
            remove.push_back(i);
    }

    int cycleBound = 0;
    std::vector<int> visited(G.size());
    for (int i = 0; i < n; i++) {
        if (G.valid_vertex(i) && !visited[i] && cover[i] + revcover[i] == 1) {
            assert(matching[i] != INT_MAX && revmatching[i] != INT_MAX);
            int length = 0;
            int x = i;
            do {
                visited[x] = 1;
                x = matching[x];
                length++;
            } while (x != i);
            cycleBound += length / 2;
        }
    }
    return cycleBound;
}

} // namespace algorithms
