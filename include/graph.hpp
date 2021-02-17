#pragma once

#include <bitset>
#include <iostream>
#include <numeric>
#include <stack>
#include <vector>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

using boost::adaptors::filtered;
using boost::adaptors::transformed;

using vertex_type = int;

class Graph {
    int vertex_count;
    int added_count;
    int removed_count;
    std::stack<std::vector<int>> unified_vertices;

    std::vector<int> vertices;
    std::vector<std::vector<int>> G;
    std::vector<bool> removed;

public:
    Graph(const std::vector<std::vector<int>> &g);
    Graph(const Graph &x) = default;
    Graph(Graph &&x) = default;
    Graph() = delete;

    void add_edge(int x, int y);

    bool is_removed(int x) const;

    bool has_neighbor(int x, int y) const;
    int count_neighbors(int x) const;
    const std::vector<int> &raw_neighbors(int x) const { return G[x]; }
    auto get_neighbors(int x) const {
        return G[x] | filtered([&](const int &edge) { return !(removed[edge]); });
    }
    auto get_vertices() const { return vertices; }
    auto get_valid_vertices() const {
        return vertices | filtered([&](const int &vertex) { return valid_vertex(vertex); });
    }
    void remove_vertex(int x);
    void unremove_vertex(int x);
    bool valid_vertex(int x) const;
    int size() const;
    int num_of_vertices() const;
    bool empty_graph() const;
    void print_graph(std::ostream &out) const;
    void ununify_last();
    int unify_vertices(const int &x, const int &y, const int &v);
    void remove_edge(const int x, const int y);

    std::vector<std::vector<int>> get_g() { return G; }
};

std::ostream &operator<<(std::ostream &out, const Graph &g);
