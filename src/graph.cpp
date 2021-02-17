#include <algorithm>
#include <cassert>
#include <functional>
#include <graph.hpp>
#include <iostream>

void Graph::add_edge(int x, int y) {
    assert(!is_removed(x));
    assert(!is_removed(y));

    auto xit = std::lower_bound(G[x].begin(), G[x].end(), y);
    if (xit != G[x].end())
        assert(*xit != y);
    G[x].insert(xit, y);

    auto yit = std::lower_bound(G[y].begin(), G[y].end(), x);
    if (yit != G[y].end())
        assert(*yit != x);
    G[y].insert(yit, x);
}

int Graph::count_neighbors(int x) const {
    assert(valid_vertex(x));

    return count_if(G[x].begin(), G[x].end(), [&](const int &edge) { return !(is_removed(edge)); });
}

bool Graph::is_removed(int x) const {
    assert(x < static_cast<int>(removed.size()));
    return removed[x];
}

bool Graph::has_neighbor(int x, int y) const {
    if (is_removed(x) || is_removed(y))
        return false;

    auto it = lower_bound(G[x].begin(), G[x].end(), y);

    if (it == G[x].end() || (*it) != y || is_removed(*it))
        return false;

    return true;
}

Graph::Graph(const std::vector<std::vector<int>> &g)
    : vertex_count(g.size()), added_count(0), removed_count(0), vertices(g.size()), G(g.size()),
      removed(g.size()) {
    std::iota(vertices.begin(), vertices.end(), 0);

    for (int i = 0; i < vertex_count; i++) {
        for (const auto &j : g[i]) {
            G[i].push_back(j);
            G[j].push_back(i);
        }
    }

    for (int i = 0; i < vertex_count; i++) {
        sort(G[i].begin(), G[i].end());
        G[i].erase(unique(G[i].begin(), G[i].end()), G[i].end());
    }
}

void Graph::remove_vertex(int x) {
    if (x >= vertex_count + added_count)
        throw std::invalid_argument("Vertex indicator out of range " + std::to_string(x));

    if (is_removed(x))
        throw std::invalid_argument("Vertex already removed " + std::to_string(x));

    removed[x] = 1;
    removed_count++;
}

void Graph::unremove_vertex(int x) {
    if (x >= vertex_count + added_count)
        throw std::invalid_argument("Vertex indicator out of range " + std::to_string(x));

    if (!is_removed(x))
        throw std::invalid_argument("Vertex " + std::to_string(x) + " not removed");

    removed[x] = 0;
    removed_count--;
}

bool Graph::valid_vertex(int x) const {
    if (x >= vertex_count + added_count || is_removed(x))
        return false;

    return true;
}

int Graph::size() const {
    return vertex_count + added_count;
}

int Graph::num_of_vertices() const {
    return vertex_count + added_count - removed_count;
}

bool Graph::empty_graph() const {
    return num_of_vertices() == 0;
}

std::ostream &operator<<(std::ostream &out, const Graph &g) {
    g.print_graph(out);
    return out;
}

void Graph::print_graph(std::ostream &out) const {
    out << "Size " << vertex_count << ", added " << added_count << ", removed " << removed_count
        << "\n";

    for (int x = 0; x < vertex_count + added_count; x++) {
        if (valid_vertex(x)) {
            out << x << ": ";

            auto v = get_neighbors(x);

            for (auto nx : v)
                out << nx << " ";

            out << std::endl;
        }
    }

    out << std::endl;
}

int Graph::unify_vertices(const int &x, const int &y, const int &v) {
    assert(valid_vertex(x));
    assert(valid_vertex(y));
    assert(valid_vertex(v));

    int new_vertex = vertex_count + added_count++;
    G.push_back(std::vector<int>());
    vertices.push_back(vertices.back() + 1);
    removed.push_back(false);
    std::vector<int> &adjacency_new_vertex = G.back();

    const auto &adjacency_x = G[x];
    const auto &adjacency_y = G[y];

    int it1 = 0, it2 = 0;

    while (it1 < static_cast<int>(adjacency_x.size())
           || it2 < static_cast<int>(adjacency_y.size())) {
        if (it2 == static_cast<int>(adjacency_y.size())
            || (it1 < static_cast<int>(adjacency_x.size())
                && adjacency_x[it1] < adjacency_y[it2])) {
            if (adjacency_x[it1] != y && adjacency_x[it1] != v) {
                adjacency_new_vertex.push_back(adjacency_x[it1]);
                G[adjacency_x[it1]].push_back(new_vertex);
            }

            it1++;
        } else if (it1 == static_cast<int>(adjacency_x.size())
                   || adjacency_y[it2] < adjacency_x[it1]) {
            if (adjacency_y[it2] != x && adjacency_y[it2] != v) {
                adjacency_new_vertex.push_back(adjacency_y[it2]);
                G[adjacency_y[it2]].push_back(new_vertex);
            }

            it2++;
        } else {
            if (adjacency_x[it1] != v) {
                adjacency_new_vertex.push_back(adjacency_x[it1]);
                G[adjacency_x[it1]].push_back(new_vertex);
            }
            it1++;
            it2++;
        }
    }

    remove_vertex(x);
    remove_vertex(y);
    remove_vertex(v);
    unified_vertices.push(std::vector<int>{x, y, v});

    return new_vertex;
}

void Graph::ununify_last() // assert stack not empty..
{
    for (const auto &x : G.back())
        G[x].pop_back();

    G.pop_back();
    vertices.pop_back();
    removed.pop_back();
    added_count--;

    unremove_vertex(unified_vertices.top()[0]);
    unremove_vertex(unified_vertices.top()[1]);
    unremove_vertex(unified_vertices.top()[2]);

    unified_vertices.pop();
}

void Graph::remove_edge(const int x, const int y) {
    assert(!is_removed(x));
    assert(!is_removed(y));

    auto xit = std::lower_bound(G[x].begin(), G[x].end(), y);
    assert(xit != G[x].end());
    assert(*xit == y);
    if (std::next(xit) != G[x].end())
        assert(*std::next(xit) != y);
    G[x].erase(xit);

    auto yit = std::lower_bound(G[y].begin(), G[y].end(), x);
    assert(yit != G[y].end());
    assert(*yit == x);
    if (std::next(yit) != G[y].end())
        assert(*std::next(yit) != x);
    G[y].erase(yit);
}
