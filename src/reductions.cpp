#include <algorithms.hpp>
#include <bnb.hpp>

namespace vc_bnb {

void Solver::reduction_degree_zero() {
    std::vector<int> vertices_to_remove;
    algorithms::get_isolated_vertices(graph_, vertices_to_remove);

    if (vertices_to_remove.size() > 0)
        for (auto vertex_to_remove : vertices_to_remove) {
            if (!graph_.valid_vertex(vertex_to_remove))
                continue;

            unpick_vertex(vertex_to_remove);
        }
}

void Solver::reduction_degree_one() {
    while (true) {
        int vertex_count = graph_.num_of_vertices();
        auto adjacency_v = graph_.get_valid_vertices();
        for (auto i : adjacency_v) {
            if (graph_.count_neighbors(i) == 1)
                unpick_vertex(i);
        }

        if (vertex_count == graph_.num_of_vertices())
            return;
    }
}

void Solver::reduction_degree_two() {
    while (true) {
        int vertex_count = graph_.num_of_vertices();

        std::vector<int> vertices = graph_.get_vertices();
        for (auto vertex : vertices) {
            if (!graph_.valid_vertex(vertex) || graph_.count_neighbors(vertex) != 2)
                continue;

            auto neighbor = graph_.get_neighbors(vertex).begin();
            int u = *neighbor++;
            int w = *neighbor;

            assert(graph_.valid_vertex(u));
            assert(graph_.valid_vertex(w));

            if (states_.get_state(u) != vc_bnb::action_type::invalid
                || states_.get_state(w) != vc_bnb::action_type::invalid
                || states_.get_state(vertex) != vc_bnb::action_type::invalid
                || graph_.has_neighbor(u, w))
                continue;

            if (flags.packing) {
                disable_constraints(vertex);
                disable_constraints(u);
                disable_constraints(w);
            }

            auto unified_vertex = graph_.unify_vertices(u, w, vertex);

            if (flags.verbose)
                std::cerr << "Deg2 reduction: removing v: " << vertex << ", unifying u:" << u
                          << " & w: " << w << ", adding: " << unified_vertex << std::endl;

            if (solver_constructed_at != -1) {
                assert(static_cast<int>(x.size()) == unified_vertex);
                assert(x[vertex]);
                assert(x[u]);
                assert(x[w]);
                x.push_back(-x[vertex]);

                // We encode the fact that either
                // * v is taken, but u and w are not taken, or
                // * v is not taken, but u and w are both taken.
                // In other words, x[v] is equivalent to -x[u] and also to -x[w].
                // Since the "positive" clauses already exist (due to edges),
                // we only need to add negative clauses.
                if (solver_constructed_at == static_cast<int>(decisions_.size())) {
                    sink.emit(std::array<atom, 2>{-x[vertex], -x[u]});
                    sink.emit(std::array<atom, 2>{-x[vertex], -x[w]});
                } else {
                    auto disabler = sink.allocate();
                    sink.emit(std::array<atom, 3>{-disabler, -x[vertex], -x[u]});
                    sink.emit(std::array<atom, 3>{-disabler, -x[vertex], -x[w]});
                    new_folds.push_back(disabler);
                }
            }

            states_.add(unified_vertex);
            if (flags.packing)
                vertex_to_constraints.push_back(std::vector<int>());
            got_unified_.push_back(false);

            unified_vertices_.push_back(std::make_tuple(u, w, vertex));

            action act(unified_vertex, vc_bnb::action_type::unify);
            act.central = vertex;
            act.contracted[0] = u;
            act.contracted[1] = w;
            trail_.push(act);
            trail_.push(action(u, vc_bnb::action_type::unified));
            trail_.push(action(vertex, vc_bnb::action_type::unified));
            trail_.push(action(w, vc_bnb::action_type::unified));

            states_.update_state(u, vc_bnb::action_type::unified);
            states_.update_state(vertex, vc_bnb::action_type::unified);
            states_.update_state(w, vc_bnb::action_type::unified);
            ++current_cardinality;

            got_unified_[u] = true;
            got_unified_[vertex] = true;
            got_unified_[w] = true;
        }

        if (vertex_count == graph_.num_of_vertices())
            return;
    }
}

void Solver::reduction_domination() {
    std::vector<bool> marked_vertices(graph_.size(), false);

    std::vector<vertex_type> vertices;
    for (auto v : graph_.get_valid_vertices())
        vertices.push_back(v);
    for (auto v : vertices) {
        if (!graph_.valid_vertex(v))
            continue;

        std::fill(marked_vertices.begin(), marked_vertices.end(), false);
        for (auto u : graph_.get_neighbors(v))
            marked_vertices[u] = true;

        for (auto u : graph_.get_neighbors(v)) {
            bool dominating = true;
            for (auto w : graph_.get_neighbors(u)) {
                if (w != v && !marked_vertices[w]) {
                    dominating = false;
                    break;
                }
            }

            if (dominating) {
                if (flags.verbose)
                    std::cerr << "Dom: ";
                pick_vertex(v);
                break;
            }
        }
    }
}

void Solver::reduction_simplicial() {
    std::vector<char> marked;
    std::vector<vertex_type> neighbors;
    marked.resize(graph_.size(), 0);

    auto check_simplicial = [&](vertex_type v) -> bool {
        for (auto w : neighbors) {
            int nx = 0;
            for (auto x : graph_.get_neighbors(w)) {
                if (marked[x])
                    ++nx;
            }
            if (nx != static_cast<int>(neighbors.size()) - 1)
                return false;
        }
        return true;
    };

    std::vector<vertex_type> vertices;
    for (auto v : graph_.get_valid_vertices())
        vertices.push_back(v);
    for (auto v : vertices) {
        if (!graph_.valid_vertex(v))
            continue;

        // Clear state.
        for (auto v : neighbors)
            marked[v] = 0;
        neighbors.clear();

        for (auto w : graph_.get_neighbors(v)) {
            marked[w] = 1;
            neighbors.push_back(w);
        }

        if (check_simplicial(v)) {
            if (flags.verbose)
                std::cerr << "Simplicial: ";
            unpick_vertex(v);
        }
    }
}

void Solver::reduction_unconfined() {
    std::vector<vertex_type> stack;
    std::vector<char> s_mark;
    std::vector<char> ns_mark;
    s_mark.resize(graph_.size());
    ns_mark.resize(graph_.size());

    auto check_u = [&](vertex_type u) -> bool {
        // There must be exactly one neighbor in S.
        // There must be at most one neighbor not in N[S].
        bool seen_s = false;
        bool not_seen_ns = false;
        for (auto x : graph_.get_neighbors(u)) {
            if (s_mark[x]) {
                if (seen_s)
                    return false;
                seen_s = true;
            }
            if (!ns_mark[x]) {
                if (not_seen_ns)
                    return false;
                not_seen_ns = true;
            }
        }
        assert(seen_s);

        return true;
    };

    auto find_w = [&](vertex_type u) -> vertex_type {
        vertex_type z = -1;
        for (auto x : graph_.get_neighbors(u)) {
            if (ns_mark[x])
                continue;
            assert(!s_mark[x]);
            assert(z < 0);
            z = x;
        }

        return z;
    };

    auto is_unconfined = [&](vertex_type v) -> bool {
        while (true) {
            s_mark[v] = 1;
            ns_mark[v] = 1;
            stack.push_back(v);
            for (auto x : graph_.get_neighbors(v)) {
                if (ns_mark[x])
                    continue;
                ns_mark[x] = 1;
                stack.push_back(x);
            }

            vertex_type w = -1;
            for (auto u : stack) {
                if (s_mark[u])
                    continue;
                assert(ns_mark[u]);

                if (!check_u(u))
                    continue;
                if (w < 0) {
                    auto z = find_w(u);
                    if (z < 0)
                        return true;
                    w = z;
                }
            }

            if (w < 0)
                return false;
            assert(!s_mark[w]);
            assert(!ns_mark[w]);

            s_mark[w] = 1;
            ns_mark[w] = 1;
            stack.push_back(w);
            for (auto x : graph_.get_neighbors(w)) {
                if (ns_mark[x])
                    continue;
                ns_mark[x] = 1;
                stack.push_back(x);
            }
        }
    };

    std::vector<vertex_type> vertices;
    for (auto v : graph_.get_valid_vertices())
        vertices.push_back(v);
    for (auto v : vertices) {
        for (auto x : stack) {
            s_mark[x] = 0;
            ns_mark[x] = 0;
        }
        stack.clear();

        if (is_unconfined(v)) {
            if (flags.verbose)
                std::cerr << "Unconfined: ";
            assert(pick_vertex(v));
        }
    }
}

void Solver::reduction_funnel() {
    std::vector<char> marked_v;
    std::vector<char> marked_u;
    std::vector<vertex_type> v_neighbors;
    std::vector<vertex_type> u_neighbors;
    marked_v.resize(graph_.size(), 0);
    marked_u.resize(graph_.size(), 0);

    std::vector<vertex_type> vertices;
    for (auto v : graph_.get_valid_vertices())
        vertices.push_back(v);
    for (auto v : vertices) {
        auto check_funnel = [&](vertex_type z) -> bool {
            for (auto w : graph_.get_neighbors(v)) {
                if (w == z)
                    continue;
                int nx = 0;
                for (auto x : graph_.get_neighbors(w)) {
                    if (!marked_v[x] || x == z)
                        continue;
                    ++nx;
                }
                // w will see neither u nor itself, hence -2.
                if (nx != static_cast<int>(v_neighbors.size()) - 2)
                    return false;
            }
            return true;
        };

        if (!graph_.valid_vertex(v))
            continue;

        // Reset state.
        for (auto w : v_neighbors)
            marked_v[w] = 0;
        v_neighbors.clear();
        for (auto w : u_neighbors)
            marked_u[w] = 0;
        u_neighbors.clear();

        // Discover neighbors of v.
        for (auto w : graph_.get_neighbors(v)) {
            marked_v[w] = 1;
            v_neighbors.push_back(w);
        }

        // We are looking for a *non-empty* clique in N(v) - {u}.
        // We could special case deg-1 here but that will produce misleading statistics
        // and changes in the graph would cause deg-1 to be re-run anyway.
        if (v_neighbors.size() == 1)
            continue;

        // Try to find a u that realizes a funnel.
        vertex_type u = -1;
        for (auto z : graph_.get_neighbors(v)) {
            if (!check_funnel(z))
                continue;
            u = z;
            break;
        }
        if (u < 0)
            continue;

        // Discover neighbors of u.
        for (auto w : graph_.get_neighbors(u)) {
            marked_u[w] = 1;
            u_neighbors.push_back(w);
        }

        // Remove the funnel from the graph.
        action funnel_act(action_type::funnel);
        funnel_act.contracted[0] = v;
        funnel_act.contracted[1] = u;
        trail_.push(funnel_act);
        states_.update_state(v, action_type::funnel);
        states_.update_state(u, action_type::funnel);

        disable_constraints(v);
        disable_constraints(u);
        graph_.remove_vertex(v);
        graph_.remove_vertex(u);
        ++current_cardinality;

        if (solver_constructed_at != -1) {
            assert(x[v]);
            assert(x[u]);
            if (solver_constructed_at == static_cast<int>(decisions_.size())) {
                sink.emit(std::array<atom, 2>{-x[v], -x[u]});
            } else {
                auto disabler = sink.allocate();
                sink.emit(std::array<atom, 3>{-disabler, -x[v], -x[u]});
                new_folds.push_back(disabler);
            }
        }

        // Pick common neighbors.
        for (vertex_type w : v_neighbors) {
            if (!marked_u[w])
                continue;
            assert(w != v && w != u);
            if (flags.verbose)
                std::cerr << "Funnel: ";
            pick_vertex(w);
        }
        for (vertex_type w : u_neighbors) {
            if (!marked_v[w])
                continue;
            assert(!graph_.valid_vertex(w));
        }

        for (auto x : v_neighbors) {
            if (marked_u[x] || x == u)
                continue;
            for (auto y : u_neighbors) {
                if (marked_v[y] || y == v)
                    continue;
                if (graph_.has_neighbor(x, y))
                    continue;

                action edge_act{action_type::add_edge};
                edge_act.contracted[0] = x;
                edge_act.contracted[1] = y;
                trail_.push(edge_act);
                graph_.add_edge(x, y);
                assert(graph_.has_neighbor(x, y));
            }
        }
    }
}

bool Solver::clean_up() {
    bool close_branch = false;

    bool changes_zero = false;
    bool changes_one = false;
    bool changes_two = false;
    bool changes_dom = false;
    bool changes_simplicial = false;
    bool changes_unconfined = false;
    bool changes_funnel = false;

    while (true) {
        int vertex_count = graph_.num_of_vertices();
        int total_vertex_count = graph_.size();

        auto start = std::chrono::high_resolution_clock::now();
        reduction_degree_zero();
        auto end = std::chrono::high_resolution_clock::now();
        stats_.reduction_zero_count += vertex_count - graph_.num_of_vertices();
        stats_.reduction_zero_time += end - start;
        if (vertex_count - graph_.num_of_vertices() > 0)
            changes_zero = true;

        if (flags.degree1) {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_degree_one();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_one_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_one_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_one = true;
        }

        if (flags.domination) {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_domination();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_dom_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_dom_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_dom = true;
        }

        {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_simplicial();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_simplicial_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_simplicial_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_simplicial = true;
        }

        {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_unconfined();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_unconfined_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_unconfined_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_unconfined = true;
        }

        if (decisions_.size() == 1 || flags.use_funnel) {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_funnel();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_funnel_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_funnel_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_funnel = true;
        }

        if (flags.deg2) {
            int vertex_count_ = graph_.num_of_vertices();
            auto start = std::chrono::high_resolution_clock::now();
            reduction_degree_two();
            auto end = std::chrono::high_resolution_clock::now();
            stats_.reduction_two_count += vertex_count_ - graph_.num_of_vertices();
            stats_.reduction_two_time += end - start;
            if (vertex_count_ - graph_.num_of_vertices() > 0)
                changes_two = true;
        }

        if (flags.packing) {
            if (!check_constraints()) {
                close_branch = true;
                break;
            }
        }

        if (total_vertex_count == graph_.size() && vertex_count == graph_.num_of_vertices())
            break;
    }

    if (changes_zero)
        stats_.reduction_zero_branch_count++;
    if (changes_one)
        stats_.reduction_one_branch_count++;
    if (changes_two)
        stats_.reduction_two_branch_count++;
    if (changes_dom)
        stats_.reduction_dom_branch_count++;
    if (changes_simplicial)
        stats_.reduction_simplicial_branch_count++;
    if (changes_unconfined)
        stats_.reduction_unconfined_branch_count++;
    if (changes_funnel)
        stats_.reduction_funnel_branch_count++;

    return close_branch;
}

int Solver::lp_reduction(int lowerbound) {
    std::vector<int> add_vertices;
    std::vector<int> remove_vertices;
    lowerbound = std::max(
        lowerbound, algorithms::cycleCoverILPRelaxation(graph_, add_vertices, remove_vertices));

    for (int i = 0; i < static_cast<int>(add_vertices.size()); i++) {
        if (flags.verbose)
            std::cerr << "lp: ";
        pick_vertex(add_vertices[i]);
    }

    for (int i = 0; i < static_cast<int>(remove_vertices.size()); i++) {
        if (flags.verbose)
            std::cerr << "lp: ";
        unpick_vertex(remove_vertices[i]);
    }
    return lowerbound;
}

} // namespace vc_bnb
