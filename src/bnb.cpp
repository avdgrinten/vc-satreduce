#include <action.hpp>
#include <algorithm>
#include <algorithms.hpp>
#include <array>
#include <bnb.hpp>
#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <numeric>
#include <sat.hpp>
#include <signal.h>
#include <states.hpp>
#include <typeinfo>

using boost::adaptors::filtered;
using boost::adaptors::transformed;

/* flag for output stats after signal */
volatile sig_atomic_t timeout_signal = 0;

static void catch_sigxcpu(int sig_nr) {
    timeout_signal = 1;
}

static constexpr bool encode_immediately = true;
static constexpr bool use_probing = false;
static constexpr bool do_exclusion_first = true;
static constexpr bool use_cone_of_influence = true;

namespace vc_bnb {
/*
 * Implementation of solver
 */

Solver::Solver(Graph &graph)
    : best_solution(graph.size()), graph_(graph), trail_(), start_vertex_count(graph.size()),
      states_(graph.size()), vertex_to_constraints(graph.size()), got_unified_(graph.size()) {}

bool Solver::is_picked(int vertex) {
    return states_.get_state(vertex) == action_type::pick;
}

bool Solver::is_unpicked(int vertex) {
    return states_.get_state(vertex) == action_type::unpick;
}

bool Solver::is_removed(int vertex) const {
    return !graph_.valid_vertex(vertex);
}

bool Solver::is_unified(int vertex) const {
    return vertex >= start_vertex_count;
}

bool Solver::got_unified(int vertex) const {
    assert(vertex < static_cast<int>(got_unified_.size()));
    return got_unified_[vertex];
}

bool Solver::is_constraint_removed(int constraint_id) const {
    assert(constraint_id < static_cast<int>(constraint_removed.size()));
    return constraint_removed[constraint_id];
}

std::tuple<int, int, int> Solver::get_unified_vertices(int vertex) const {
    assert(vertex - start_vertex_count < static_cast<int>(unified_vertices_.size()));
    return unified_vertices_[vertex - start_vertex_count];
}

void Solver::project_solution(Solution &out, const std::vector<vertex_type> &vc) {
    enum class state { invalid, included, excluded };

    std::vector<state> s;
    s.resize(graph_.size(), state::invalid);
    for (vertex_type x : graph_.get_valid_vertices())
        s[x] = state::excluded;

    for (vertex_type x : vc) {
        assert(s[x] == state::excluded);
        assert(graph_.valid_vertex(x));
        s[x] = state::included;
    }

    for (int i = trail_.length() - 1; i >= 0; --i) {
        auto act = trail_.action_at(i);
        if (act.type == action_type::pick) {
            auto x = act.vertex;
            assert(s[x] == state::invalid);
            s[x] = state::included;
        } else if (act.type == action_type::unpick) {
            auto x = act.vertex;
            assert(s[x] == state::invalid);
            s[x] = state::excluded;
        } else if (act.type == action_type::unify) {
            auto u = act.vertex;
            auto x = act.central;
            auto v = act.contracted[0];
            auto w = act.contracted[1];
            assert(s[x] == state::invalid);
            assert(s[v] == state::invalid);
            assert(s[w] == state::invalid);
            if (s[u] == state::included) {
                s[x] = state::excluded;
                s[v] = state::included;
                s[w] = state::included;
            } else {
                assert(s[u] == state::excluded);
                s[x] = state::included;
                s[v] = state::excluded;
                s[w] = state::excluded;
            }
            s[u] = state::invalid;
        } else if (act.type == action_type::funnel) {
            auto v = act.contracted[0];
            auto u = act.contracted[1];
            assert(s[v] == state::invalid);
            assert(s[u] == state::invalid);
            bool all_v_neighbors = true;
            bool all_u_neighbors = true;
            for (auto w : graph_.raw_neighbors(v)) {
                if (s[w] == state::invalid)
                    continue;
                if (s[w] != state::included)
                    all_v_neighbors = false;
            }
            for (auto w : graph_.raw_neighbors(u)) {
                if (s[w] == state::invalid)
                    continue;
                if (s[w] != state::included)
                    all_u_neighbors = false;
            }
            if (all_v_neighbors) {
                s[u] = state::included;
                s[v] = state::excluded;
            } else {
                assert(all_u_neighbors);
                s[v] = state::included;
                s[u] = state::excluded;
            }
        } else {
            assert(act.type == action_type::add_edge || act.type == action_type::unified
                   || act.type == action_type::constraint_created
                   || act.type == action_type::constraint_removed
                   || act.type == action_type::removed_from_constraint_pick
                   || act.type == action_type::removed_from_constraint_unpick
                   || act.type == action_type::constraint_disabled);
        }
    }

    // TODO: check that we only include original vertices.
    out = Solution{graph_.size()};
    for (int x = 0; x < graph_.size(); ++x) {
        if (s[x] == state::included)
            out.select_vertex(x);
    }
}

int Solver::get_next_vertex() {
    int next_vertex = INT_MAX;
    int next_vertex_neighbors_count = -1;

    for (auto vertex : graph_.get_valid_vertices()) {
        int vertex_neighbors_count = graph_.count_neighbors(vertex);
        if (vertex_neighbors_count > next_vertex_neighbors_count) {
            next_vertex = vertex;
            next_vertex_neighbors_count = vertex_neighbors_count;
        }
    }

    return next_vertex;
}

void Solver::branching() {
    auto next_vertex = get_next_vertex();
    assert(graph_.valid_vertex(next_vertex));

    found_mirrors.clear();
    found_satellites.clear();

    auto start = std::chrono::high_resolution_clock::now();
    if (!flags.packing && flags.mirror_branch && flags.satellite_branch) {
        algorithms::find_mirrors_and_satellites(graph_, next_vertex, found_mirrors,
                                                found_satellites);
    } else if (flags.mirror_branch) {
        algorithms::find_mirrors(graph_, next_vertex, found_mirrors);
    } else if (flags.satellite_branch) {
        algorithms::find_satellites(graph_, next_vertex, found_satellites);
    }
    auto end = std::chrono::high_resolution_clock::now();
    stats_.satmir_time += end - start;

    auto current_trail = trail_.length();
    auto inherited_bound = decisions_.back().bound;

    decisions_.emplace_back();
    auto &d = decisions_.back();
    d.vertex = next_vertex;
    d.trail_ptr = current_trail;
    d.bound = inherited_bound;

    // Push a stack canary for sanity checking.
    alternatives.push_back(-1);

    if (!found_mirrors.empty()) {
        d.alternative = alternative_type::exclude;

        // Pick the vertex and all of its mirrors.
        assert(graph_.valid_vertex(next_vertex));
        pick_vertex(next_vertex);
        for (vertex_type v : found_mirrors)
            pick_vertex(v);
    } else {
        if (do_exclusion_first) {
            d.alternative = alternative_type::include;

            if (flags.packing && flags.unpick_constraint)
                create_unpick_constraint(next_vertex);

            // Exclude the vertex and all satellites.
            assert(graph_.valid_vertex(next_vertex));
            unpick_vertex(next_vertex);
            for (vertex_type u : found_satellites)
                unpick_vertex(u);
        } else {
            d.alternative = alternative_type::exclude_satellites;
            for (vertex_type u : found_satellites)
                alternatives.push_back(u);

            if (flags.packing && flags.pick_constraint)
                create_pick_constraint(next_vertex);

            assert(graph_.valid_vertex(next_vertex));
            pick_vertex(next_vertex);
        }
    }

    if (!found_mirrors.empty()) {
        stats_.mirror_count += found_mirrors.size();
        ++stats_.mirror_branches;
    }
    if (!found_satellites.empty()) {
        stats_.satellite_count += found_satellites.size();
        ++stats_.satellite_branches;
    }
}

void Solver::track_decision(action act) {
    states_.update_state(act.vertex, act.type);
    trail_.push(act);
}

bool Solver::pick_vertex(int vertex) {
    if (!graph_.valid_vertex(vertex)) {
        assert(is_picked(vertex));
        return false;
    }

    track_decision(action{vertex, action_type::pick});

    if (flags.verbose) {
        std::cerr << "Picking vertex " << vertex << std::endl;
    }

    if (solver_constructed_at != -1) {
        if (solver_constructed_at == static_cast<int>(decisions_.size())) {
            assert(x[vertex]);
            sink.emit(std::array<atom, 1>{x[vertex]});
        } else {
            new_picked_vertices.push_back(vertex);
        }
    }

    graph_.remove_vertex(vertex);
    ++current_cardinality;
    ++changed_vertices;

    if (flags.packing) {
        if (!remove_in_constraints(vertex, true))
            return false;
    }

    return true;
}

bool Solver::unpick_vertex(int vertex) {
    if (!graph_.valid_vertex(vertex)) {
        assert(is_unpicked(vertex));
        return false;
    }

    track_decision(action{vertex, action_type::unpick});

    if (flags.verbose) {
        std::cerr << "Unpicking vertex " << vertex << std::endl;
    }

    if (solver_constructed_at != -1) {
        if (solver_constructed_at == static_cast<int>(decisions_.size())) {
            assert(x[vertex]);
            sink.emit(std::array<atom, 1>{-x[vertex]});
        } else {
            new_unpicked_vertices.push_back(vertex);
        }
    }

    graph_.remove_vertex(vertex);
    ++changed_vertices;

    if (flags.packing)
        if (!remove_in_constraints(vertex, false))
            return false;

    // pick neighbors
    auto adjacency_v = graph_.get_neighbors(vertex);
    for (auto neighbor : adjacency_v) {
        if (flags.verbose)
            std::cerr << "neighbor ";

        pick_vertex(neighbor);
    }

    return true;
}

bool Solver::unpick_satellites() {
    auto satellites = satellites_.back();
    satellites_.pop_back();

    for (auto satellite_ : satellites) {
        if (!graph_.valid_vertex(satellite_)
            || states_.get_state(satellite_) == action_type::unpick)
            continue;

        if (flags.verbose)
            std::cerr << "unpick satellite: " << satellite_ << std::endl;

        unpick_vertex(satellite_);
    }
    return true;
}

// removes variable of vertex from left side of any constraint it occurs in
bool Solver::remove_in_constraints(int vertex, bool decrease_right_side) {
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < vertex_to_constraints[vertex].size(); i++) {
        auto constraint_id = vertex_to_constraints[vertex][i];
        assert(constraint_id < static_cast<int>(constraints.size()));

        if (constraint_removed[constraint_id] || !constraint_enabled[constraint_id])
            continue;

        action_type type;
        if (decrease_right_side)
            type = action_type::removed_from_constraint_pick;
        else
            type = action_type::removed_from_constraint_unpick;
        trail_.push(action(vertex, type, constraint_id));

        auto &constraint = constraints[constraint_id];

        constraint.update(vertex, decrease_right_side);
    }

    auto end = std::chrono::high_resolution_clock::now();
    stats_.packing_time += end - start;
    return true;
}

bool Solver::check_constraints() {
    std::vector<char> marked;
    std::vector<int> neighbor_cnt;
    std::vector<vertex_type> visted;
    marked.resize(graph_.size(), 0);
    neighbor_cnt.resize(graph_.size(), 0);

    auto reduce_zero_rhs = [&](size_t cidx) -> bool {
        // Clean the state.
        for (vertex_type x : visted) {
            marked[x] = false;
            neighbor_cnt[x] = 0;
        }
        visted.clear();

        auto &constraint = constraints[cidx];

        for (auto v : constraint.get_left_side()) {
            marked[v] = true;
            visted.push_back(v);
        }

        for (auto v : constraint.get_left_side()) {
            for (auto u : graph_.get_neighbors(v)) {
                if (marked[u])
                    continue;
                neighbor_cnt[u]++;
                if (neighbor_cnt[u] == 1)
                    visted.push_back(u);
            }
        }

        for (auto u : visted) {
            if (marked[u])
                continue;
            assert(neighbor_cnt[u]);

            if (neighbor_cnt[u] == 1) {
                stats_.packing_reduction2_count++;
                assert(graph_.valid_vertex(u));

                // Create new constraint \sum N^+(u) ≤ |N^+(u)| - 1.
                std::vector<int> adj_vector;
                for (int w : graph_.get_neighbors(u)) {
                    if (marked[w] || neighbor_cnt[w])
                        continue;
                    adj_vector.push_back(w);
                }

                auto rhs = adj_vector.size() - 1;
                create_constraint(std::move(adj_vector), rhs);
            }
        }

        {
            // Need to get the constraint again, since create_constraint() above
            // resizes the constraint vector.
            auto &constraint = constraints[cidx];

            if (constraint.get_left_side().size() == 1) {
                stats_.packing_reduction1_count++;
                unpick_vertex(constraint.get_left_side().back());
            } else if (constraint.get_left_side().size() > 1) {
                const auto &lhs = constraint.get_left_side();
                for (int i = 0; i < static_cast<int>(lhs.size()) - 1; i++) {
                    for (int j = i + 1; j < static_cast<int>(lhs.size()); j++) {
                        if (graph_.has_neighbor(lhs[i], lhs[j]))
                            return false;
                    }
                }

                stats_.packing_reduction2_count++;
                for (auto v : lhs)
                    unpick_vertex(v);
            }
            remove_constraint(cidx);
        }

        return true;
    };

    auto reduce_nonzero_rhs = [&](size_t cidx) {
        // Clean the state.
        for (vertex_type x : visted) {
            marked[x] = false;
            neighbor_cnt[x] = 0;
        }
        visted.clear();

        auto &constraint = constraints[cidx];

        for (auto v : constraint.get_left_side()) {
            marked[v] = true;
            visted.push_back(v);
        }

        for (auto v : constraint.get_left_side()) {
            for (auto u : graph_.get_neighbors(v)) {
                if (marked[u])
                    continue;
                neighbor_cnt[u]++;
                if (neighbor_cnt[u] == 1)
                    visted.push_back(u);
            }
        }

        for (auto u : visted) {
            if (marked[u])
                continue;
            assert(neighbor_cnt[u]);

            // Need to get the constraint again, since create_constraint() below
            // resizes the constraint vector.
            auto &constraint = constraints[cidx];

            if (neighbor_cnt[u] > constraint.get_right_side()) {
                if (flags.verbose)
                    std::cerr << "PR#3: pick " << u << std::endl;
                stats_.packing_reduction3_count++;
                assert(graph_.valid_vertex(u));
                pick_vertex(u);

                // Create new constraint \sum N(u) ≤ |N(u)| - 2.
                std::vector<int> adj_vector;
                auto adj_range = graph_.get_neighbors(u);
                std::copy(adj_range.begin(), adj_range.end(), std::back_inserter(adj_vector));

                auto rhs = adj_vector.size() - 2;
                create_constraint(std::move(adj_vector), rhs);
            }
        }
    };

    for (size_t constraint_id = 0; constraint_id < constraints.size(); constraint_id++) {
        if (!constraint_enabled[constraint_id] || constraint_removed[constraint_id])
            continue;

        auto &constraint = constraints[constraint_id];

        for (vertex_type v : constraint.get_left_side()) {
            assert(graph_.valid_vertex(v));
        }

        // Close branches due to trivally unsatisfiable constraints.
        if (constraint.get_right_side() < 0) {
            if (flags.verbose)
                std::cerr << "backtrack because unsatisfied constraint (PR#0)" << std::endl;
            stats_.packing_branches_closed++;
            return false;
        }

        // Discard trivally satisfied constraints.
        if (constraint.get_left_side().empty()) {
            remove_constraint(constraint_id);
            continue;
        }
        if (static_cast<int>(constraint.get_left_side().size()) <= constraint.get_right_side()) {
            remove_constraint(constraint_id);
            continue;
        }

        if (constraint.get_right_side() == 0) {
            auto start = std::chrono::high_resolution_clock::now();
            if (!reduce_zero_rhs(constraint_id))
                return false;
            stats_.packing_rhs_zero_time += std::chrono::high_resolution_clock::now() - start;
        } else if (constraint.get_left_side().size() > 0) {
            auto start = std::chrono::high_resolution_clock::now();
            reduce_nonzero_rhs(constraint_id);
            stats_.packing_rhs_nonzero_time += std::chrono::high_resolution_clock::now() - start;
        }
    }

    return true;
}

void Solver::remove_constraint(int constraint_id) {
    assert(!constraint_removed[constraint_id]);
    constraint_removed[constraint_id] = true;
    trail_.push(action(0, action_type::constraint_removed, constraint_id));
}

void Solver::create_pick_constraint(int vertex) {
    std::vector<int> adj_vector;
    auto adj_range = graph_.get_neighbors(vertex);
    std::copy(adj_range.begin(), adj_range.end(), std::back_inserter(adj_vector));

    auto rhs = adj_vector.size() - 1;
    create_constraint(std::move(adj_vector), rhs);
}

void Solver::create_unpick_constraint(int vertex) {
    std::vector<char> mark;
    mark.resize(graph_.size(), 0);
    mark[vertex] = 1;
    for (auto w : graph_.get_neighbors(vertex))
        mark[w] = 1;

    std::vector<int> neighbors_plus;
    for (auto w : graph_.get_neighbors(vertex)) {
        neighbors_plus.clear();
        for (auto x : graph_.get_neighbors(w))
            if (!mark[x])
                neighbors_plus.push_back(x);
        create_constraint(neighbors_plus, neighbors_plus.size() - 1);
    }
}

// creates a constraint like: \sum v \in left_side ≤ right_side
void Solver::create_constraint(std::vector<int> left_side, int right_side) {
    // Otherwise the constraint is trival at creation time.
    assert(static_cast<int>(left_side.size()) > right_side);

    for (vertex_type v : left_side)
        assert(graph_.valid_vertex(v));

    auto start = std::chrono::high_resolution_clock::now();

    if (left_side.size()) {
        // Get the LHS vertex that occurs least often in constraints,
        // check redundancy against all constraints that contain it.
        auto vit = std::min_element(
            left_side.begin(), left_side.end(), [&](vertex_type xl, vertex_type xr) -> bool {
                return vertex_to_constraints[xl].size() < vertex_to_constraints[xr].size();
            });
        assert(vit != left_side.end());
        for (auto cidx : vertex_to_constraints[*vit]) {
            if (constraints[cidx].equals(left_side, right_side))
                return;
        }
    }

    constraints.emplace_back(left_side, right_side);
    int index = constraints.size() - 1;
    for (auto u : constraints[index].get_left_side())
        vertex_to_constraints[u].push_back(index);
    constraint_enabled.push_back(true);
    constraint_removed.push_back(false);
    constraint_literal.push_back(0);
    trail_.push(action(0, action_type::constraint_created, index));

    if (flags.constraint_totalizer && solver_constructed_at != -1) {
        create_totalizer_from_constraint(index);
        if (solver_constructed_at == static_cast<int>(decisions_.size())) {
            assert(constraint_literal[index]);
            sink.emit(std::array<atom, 1>{constraint_literal[index]});
        } else {
            new_constraints.push_back(index);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    stats_.packing_time += end - start;
    ++stats_.packing_constraints;
}

void Solver::disable_constraints(int vertex) {
    for (size_t i = 0; i < vertex_to_constraints[vertex].size(); i++) {
        auto constraint_id = vertex_to_constraints[vertex][i];

        if (constraint_removed[constraint_id] || !constraint_enabled[constraint_id])
            continue;

        constraint_enabled[constraint_id] = false;
        trail_.push(action(0, action_type::constraint_disabled, constraint_id));
    }
}

void Solver::undo_step() {
    auto act = trail_.pop();
    int vertex = act.vertex;

    // backtrack constraints
    if (flags.packing && act.constraint_id != -1) {
        auto constraint_id = act.constraint_id;
        if (act.type == action_type::constraint_created) {
            if (flags.constraint_totalizer && solver_constructed_at != -1) {
                assert(!new_constraints.empty() && new_constraints.back() == constraint_id);
                assert(constraint_literal[constraint_id]);
                // Emit the inverse of the constraint literal to effectively delete the constraint.
                sink.emit(std::array<atom, 1>{-constraint_literal[constraint_id]});
                new_constraints.pop_back();
            }

            for (auto v : constraints[constraint_id].get_original_left_side()) {
                for (size_t i = 0; i < vertex_to_constraints[v].size(); i++) {
                    if (vertex_to_constraints[v][i] == constraint_id) {
                        vertex_to_constraints[v].erase(vertex_to_constraints[v].begin() + i);
                        break;
                    }
                }
            }

            assert(constraints.size() > 0);
            assert(constraint_enabled.size() > 0);
            assert(constraint_removed.size() > 0);

            constraints.pop_back();
            constraint_enabled.pop_back();
            constraint_removed.pop_back();
            constraint_literal.pop_back();
        } else if (act.type == action_type::constraint_removed) {
            assert(constraint_id < static_cast<int>(constraint_removed.size()));
            assert(constraint_removed[constraint_id]);
            constraint_removed[constraint_id] = false;
        } else if (act.type == action_type::constraint_disabled) {
            assert(constraint_id < static_cast<int>(constraint_enabled.size()));
            assert(!constraint_enabled[constraint_id]);
            constraint_enabled[constraint_id] = true;
        } else if (act.type == action_type::removed_from_constraint_pick
                   || act.type == action_type::removed_from_constraint_unpick) {
            assert(!got_unified(vertex));
            assert(constraint_id < static_cast<int>(constraint_enabled.size()));
            assert(constraint_id < static_cast<int>(constraint_removed.size()));
            assert(constraint_enabled[constraint_id]);
            assert(!constraint_removed[constraint_id]);

            auto &constraint = constraints[constraint_id];
            bool decreased_right_side = false;
            if (act.type == action_type::removed_from_constraint_pick)
                decreased_right_side = true;
            constraint.backtrack(vertex, decreased_right_side);
        } else {
            assert(!"Invalid action");
        }
        return;
    }

    if (flags.verbose)
        std::cerr << "Backtrack vertex: " << vertex << " action: " << act.type << std::endl;

    if (act.type == action_type::pick) {
        assert(is_removed(vertex));
        graph_.unremove_vertex(vertex);
        states_.reset_state(vertex);

        assert(current_cardinality);
        --current_cardinality;

        if (solver_constructed_at != -1) {
            assert(!new_picked_vertices.empty() && vertex == new_picked_vertices.back());
            new_picked_vertices.pop_back();
        }

        if (changed_vertices > 0)
            changed_vertices--;
    } else if (act.type == action_type::unpick) {
        assert(is_removed(vertex));
        graph_.unremove_vertex(vertex);
        states_.reset_state(vertex);

        if (solver_constructed_at != -1) {
            assert(!new_unpicked_vertices.empty() && vertex == new_unpicked_vertices.back());
            new_unpicked_vertices.pop_back();
        }

        if (changed_vertices > 0)
            changed_vertices--;
    } else if (act.type == action_type::unify) {
        graph_.ununify_last();
        states_.remove(act.vertex);
        if (flags.packing)
            vertex_to_constraints.pop_back();
        got_unified_.pop_back();
        assert(current_cardinality);
        --current_cardinality;

        unified_vertices_.pop_back();

        if (solver_constructed_at != -1) {
            assert(static_cast<int>(x.size()) == vertex + 1);
            assert(x[act.central]);
            assert(x[act.contracted[0]]);
            assert(x[act.contracted[1]]);
            assert(x[vertex] == -x[act.central]);
            x.pop_back();
            assert(!new_folds.empty());
            assert(new_folds.back());
            // Emit the inverse literal to effectively undo the fold.
            sink.emit(std::array<atom, 1>{-new_folds.back()});
            new_folds.pop_back();
        }
    } else if (act.type == action_type::unified) {
        assert(vertex < static_cast<int>(got_unified_.size()));
        states_.reset_state(vertex);
        got_unified_[vertex] = false;
    } else if (act.type == action_type::add_edge) {
        assert(graph_.has_neighbor(act.contracted[0], act.contracted[1]));
        graph_.remove_edge(act.contracted[0], act.contracted[1]);
        assert(!graph_.has_neighbor(act.contracted[0], act.contracted[1]));
    } else if (act.type == action_type::funnel) {
        assert(graph_.is_removed(act.contracted[0]));
        assert(graph_.is_removed(act.contracted[1]));
        states_.reset_state(act.contracted[0]);
        states_.reset_state(act.contracted[1]);
        graph_.unremove_vertex(act.contracted[0]);
        graph_.unremove_vertex(act.contracted[1]);
        --current_cardinality;

        if (solver_constructed_at != -1) {
            assert(!new_folds.empty());
            assert(new_folds.back());
            // Emit the inverse literal to effectively undo the funnel.
            sink.emit(std::array<atom, 1>{-new_folds.back()});
            new_folds.pop_back();
        }
    } else {
        assert(!"Invalid action");
    }
}

bool Solver::backtrack() {
    assert(!decisions_.empty());

    if (decisions_.size() == 1) {
        if (flags.verbose) {
            std::cerr << "No alternatives left" << std::endl;
            std::cerr << current_cardinality << std::endl;
            std::cerr << best_solution.size() << std::endl;
        }
        return false;
    }

    // Invalidate the SAT solver.
    if (solver_constructed_at == static_cast<int>(decisions_.size()))
        solver_constructed_at = -1;

    auto d = decisions_.back();
    decisions_.pop_back();
    assert(d.vertex != -1);

    while (trail_.length() > d.trail_ptr)
        undo_step();

    if (d.alternative == alternative_type::include) {
        if (flags.packing && flags.pick_constraint)
            create_pick_constraint(d.vertex);

        assert(graph_.valid_vertex(d.vertex));
        pick_vertex(d.vertex);
    } else if (d.alternative == alternative_type::exclude) {
        assert(graph_.valid_vertex(d.vertex));
        unpick_vertex(d.vertex);
    } else {
        assert(d.alternative == alternative_type::exclude_satellites);

        if (flags.packing && flags.unpick_constraint)
            create_unpick_constraint(d.vertex);

        // Exclude the vertex and all satellites.
        assert(graph_.valid_vertex(d.vertex));
        unpick_vertex(d.vertex);

        while (true) {
            assert(!alternatives.empty());
            if (alternatives.back() == -1)
                break;
            unpick_vertex(alternatives.back());
            alternatives.pop_back();
        }
    }

    // Remove the canary.
    assert(!alternatives.empty() && alternatives.back() == -1);
    alternatives.pop_back();

    return true;
}

bool Solver::solve_ccs() {
    assert(ccs_.size() > 1);

    // Pick one component that we continue to solve with this instance.
    size_t pc = 0;
    for (size_t c = 1; c < ccs_.size(); ++c) {
        if (ccs_[c].size() > ccs_[pc].size())
            pc = c;
    }

    ++stats_.cc_branches;

    // Solve other CCs by creating a new instance of the solver.
    std::vector<int> mapping(graph_.size(), INT_MAX);
    std::vector<int> backmapping;
    std::vector<std::vector<int>> adjlist;
    for (size_t c = 0; c < ccs_.size(); ++c) {
        if (c == pc)
            continue;

        int sz = 0;
        backmapping.resize(ccs_[c].size());
        std::fill(backmapping.begin(), backmapping.end(), INT_MAX);
        for (auto x : ccs_[c]) {
            // TODO: we disable the constraints here; this can sometimes
            // be avoided by propagating the constraints to the subsolver.
            disable_constraints(x);
            assert(mapping[x] == INT_MAX);
            mapping[x] = sz;
            backmapping[sz] = x;
            sz++;
        }
        assert(sz == static_cast<int>(ccs_[c].size()));
        adjlist.resize(sz);

        for (int i = 0; i < sz; i++) {
            adjlist[i].clear();
            for (auto x : graph_.get_neighbors(backmapping[i])) {
                x = mapping[x];
                adjlist[i].push_back(x);
                assert(x >= 0 && x < sz);
            }
        }

        Graph component_graph(adjlist);
        Solver componentsolver(component_graph);

        // flags propagation
        componentsolver.flags = flags;
        if (ccs_[c].size() < 50)
            componentsolver.flags.use_sat_solver = false;

        if (!componentsolver.solve())
            return false;

        auto best = componentsolver.get_best_solution();
        for (auto x : ccs_[c]) {
            if (best[mapping[x]]) {
                pick_vertex(x);
            } else {
                unpick_vertex(x);
            }
        }

        // adding stats to this solver
        stats_ = stats_ + componentsolver.stats_;
    }

    bool success = check_constraints();
    // TODO: see the TODO above.
    assert(success);

    return true;
}

// Ask SAT solver whether there is a VC of the current graph
// that satisfies: lb <= |VC| < ub.
bool Solver::search_oracle(oracle_strategy strategy, double limit_factor) {
    auto &d = decisions_.back();

    std::vector<atom> inner_t;
    if (solver_constructed_at < 0) {
        sink.reset();

        if (strategy == oracle_strategy::initial_ub) {
            sink.configure_sat();
        } else {
            sink.configure_unsat();
        }

        new_picked_vertices.clear();
        new_unpicked_vertices.clear();
        new_constraints.clear();
        new_folds.clear();
        std::fill(constraint_literal.begin(), constraint_literal.end(), 0);

        x.clear();
        c.clear(); // Array of all x_v variables.
        y.clear(); // Array of all clique variables
        t.clear(); // Output of the totalizer.

        x.resize(graph_.size());
        for (int i = 0; i < graph_.size(); i++) {
            if (graph_.valid_vertex(i))
                x[i] = sink.allocate();

            if (flags.vertex_encoding && graph_.valid_vertex(i))
                c.push_back(x[i]);
        }

        // Encode VC constraints.
        for (auto v : graph_.get_valid_vertices()) {
            for (auto u : graph_.get_neighbors(v)) {
                assert(graph_.valid_vertex(u));
                sink.emit(std::array<int, 2>{x[v], x[u]});
            }
        }

        if (!flags.vertex_encoding) {
            // Encode cliques
            algorithms::clique_cover(graph_, solver_cliques);
            if (solver_cliques.size() > 2000)
                return false;
            for (auto clique : solver_cliques) {
                std::vector<int> clause;
                for (auto v : clique) {
                    assert(graph_.valid_vertex(v));
                    clause.push_back(-x[v]);
                }

                y.push_back(sink.allocate());

                c.push_back(y.back());
                clause.push_back(c.back());
                sink.emit(clause);
            }
        }

        int delta;
        if (!flags.vertex_encoding) {
            delta = graph_.num_of_vertices() - static_cast<int>(c.size());
            // If the assert fails, the clique cover would improve the LB.
            assert(d.bound - current_cardinality >= delta);
        } else {
            delta = 0;
        }

        int tot_base = current_cardinality + delta;
        int tot_lb;
        int tot_ub;
        if (use_cone_of_influence) {
            tot_lb = d.bound;
            tot_ub = best_solution.size();
        } else {
            tot_lb = tot_base;
            tot_ub = tot_base + static_cast<int>(c.size());
        }
        assert(tot_lb >= tot_base);

        for (int i = 0; i < tot_ub - tot_lb; i++)
            t.push_back(sink.allocate());

        if (c.size() > 0) {
			totalizer(sink, c, t, tot_lb - tot_base, tot_ub - tot_base);
        }

        if (flags.constraint_totalizer) {
            for (size_t i = 0; i < constraints.size(); i++) {
                if (!constraint_enabled[i] || is_constraint_removed(i))
                    continue;
                create_totalizer_from_constraint(i);
                sink.emit(std::array<atom, 1>{constraint_literal[i]});
            }
        }

        solver_constructed_at = static_cast<int>(decisions_.size());
        cardinality_at_construction = current_cardinality;
        vertices_at_construction = graph_.num_of_vertices();
        id_limit_at_construction = graph_.size();
        lb_at_construction = d.bound;
        ub_at_construction = best_solution.size();
        ++stats_.solver_constructions;
    } else {
        assert(solver_constructed_at <= static_cast<int>(decisions_.size()));
    }

    int delta;
    if (!flags.vertex_encoding) {
        delta = vertices_at_construction - static_cast<int>(c.size());
        // If the assert fails, the clique cover would improve the LB.
        assert(d.bound - cardinality_at_construction >= delta);
    } else {
        delta = 0;
    }

    int tot_base = cardinality_at_construction + delta;
    int tot_lb;
    int tot_ub;
    if (use_cone_of_influence) {
        tot_lb = lb_at_construction;
        tot_ub = ub_at_construction;
    } else {
        tot_lb = tot_base;
        tot_ub = tot_base + static_cast<int>(c.size());
    }
    assert(tot_lb >= tot_base);

    int cur_lb = d.bound;
    int cur_ub = best_solution.size();
    assert(cur_lb >= tot_lb);
    assert(cur_ub <= tot_ub);

    std::vector<vertex_type> vertices_to_pick;
    std::vector<vertex_type> vertices_to_unpick;
    while (cur_lb < cur_ub && !timeout_signal) {
        // TODO: change direction
        int q;
        if (strategy == oracle_strategy::initial_ub || strategy == oracle_strategy::probe
            || flags.ub_to_lb)
            q = cur_ub - 1;
        else
            q = cur_lb;
        assert(cur_lb <= q);
        assert(q < cur_ub);

        if (flags.verbose)
            std::cerr << "invoking oracle " << cur_lb << " <= " << q << " < " << cur_ub
                      << std::endl;

        // If the upper bound improves, we can permanently fix variables.
        // In addition, this also permanently fixes the solution to the smaller
        // than the current best solution.
        for (auto b = cur_ub - 1; b < tot_ub; ++b) {
            auto ib = b - tot_lb;
            if (sink.fixed(t[ib]))
                continue;
            sink.emit(std::vector<atom>{-t[ib]});
        }

        for (auto v : new_picked_vertices) {
            assert(is_picked(v));
            assert(v < static_cast<int>(x.size()));
            assert(x[v]);
            sink.assume(x[v]);
        }

        for (auto v : new_unpicked_vertices) {
            assert(is_unpicked(v));
            assert(v < static_cast<int>(x.size()));
            assert(x[v]);
            sink.assume(-x[v]);
        }

        for (auto l : new_folds) {
            assert(l);
            sink.assume(l);
        }

        for (auto cidx : new_constraints) {
            // TODO: SAT can handle constraints that contain unified vertices just fine.
            //       We could enable them here but we would need to guarantee
            //       that they are actually generated.
            if (!constraint_enabled[cidx])
                continue;
            assert(constraint_literal[cidx]);
            sink.assume(constraint_literal[cidx]);
        }

        int iq = q - tot_lb;
        if (flags.ub_to_lb) {
            if (iq < static_cast<int>(t.size())
                && sink.fixed(t[iq]) == 0) // ignores first iteration
                sink.emit(std::vector<atom>{-t[iq]});
        } else {
            // Generates a (lhs <= m) constraint.
            auto assume_rhs = [&](int m) {
                for (int i = m; i < static_cast<int>(t.size()); i++)
                    sink.assume(-t[i]);
            };

            assume_rhs(iq);
        }

        int limit = flags.sat_conflict_limit * std::min(limit_factor, 1000.0);
        if (strategy == oracle_strategy::probe) {
            limit = 16;
        } else if (strategy == oracle_strategy::exhaustive) {
            limit = INT_MAX;
        }
        auto status = sink.solve(limit);

        if (status == sat_status::sat) // Satisfiable.
        {
            stats_.sat_count++;
            if (flags.verbose)
                std::cerr << "    VC exists" << std::endl;

            for (auto v : new_picked_vertices)
                assert(sink.value(x[v]) > 0);
            for (auto v : new_unpicked_vertices)
                assert(sink.value(x[v]) < 0);

            // Count number of totalizer variables (either clauses or vertices).
            int nc = 0;
            for (size_t i = 0; i < c.size(); ++i) {
                if (sink.value(c[i]) > 0)
                    nc++;
            }
            assert(nc <= q - tot_base);

            // Sanity check nc against the vertex cover size.
            int nv = 0;
            for (int v = 0; v < id_limit_at_construction; v++) {
                if (!x[v])
                    continue;
                if (sink.value(x[v]) > 0)
                    nv++;
            }
            assert(nv <= delta + nc); // Totalizers only guarantte <= here.

            // The cardinality enforced by the totalizer.
            int enforced_cardinality = cardinality_at_construction + delta + nc;
            assert(enforced_cardinality <= q);

            // Extract the vertex cover from the solver.
            std::vector<vertex_type> vertices;
            for (auto v : graph_.get_valid_vertices()) {
                assert(x[v]);
                if (sink.value(x[v]) > 0)
                    vertices.push_back(v);
            }

            auto cardinality = current_cardinality + static_cast<int>(vertices.size());
            if (flags.verbose)
                std::cerr << "    oracle found VC of size " << cardinality << std::endl;
            assert(cardinality <= enforced_cardinality);

            // SAT solution is always better than best solution.
            assert(cardinality < best_solution.size());
            project_solution(best_solution, vertices);
            assert(cardinality == best_solution.size());

            cur_ub = cardinality;
        } else if (status == sat_status::unsat) {
            stats_.unsat_count++;
            if (flags.verbose)
                std::cerr << "    VC does not exist" << std::endl;

            // SAT lower bound is always better than our current lower bound.
            assert(q + 1 > d.bound);
            d.bound = q + 1;

            cur_lb = q + 1;
        } else {
            assert(status == sat_status::null);
        }

        vertices_to_pick.clear();
        vertices_to_unpick.clear();
        for (auto v : graph_.get_valid_vertices()) {
            int f = sink.fixed(x[v]);
            if (f > 0) {
                vertices_to_pick.push_back(v);
                ++stats_.x_fixed;
            } else if (f < 0) {
                vertices_to_unpick.push_back(v);
                ++stats_.x_fixed;
            }
        }

        for (auto v : vertices_to_pick)
            pick_vertex(v);
        for (auto v : vertices_to_unpick) {
            // If the formula is unsatisfiable, the oracle can exclude conflicting vertices.
            if (graph_.valid_vertex(v))
                unpick_vertex(v);
        }

        if (!flags.vertex_encoding && flags.packing && flags.check_fixed_yc) {
            for (size_t i = 0; i < y.size(); i++) {
                int f = sink.fixed(y[i]);
                if (f < 0) {
                    // clique fixed with 0
                    stats_.packing_yc_fixed_count++;

                    assert(i < solver_cliques.size());
                    if (solver_cliques[i].size() == 0)
                        continue;

                    int offset = 0;
                    int unremoved_vertex = 0;
                    bool pick = false;
                    for (auto c : solver_cliques[i])
                        if (is_removed(c)) {
                            offset++;
                            if (!is_picked(c))
                                pick = true;
                        } else
                            unremoved_vertex = c;

                    if (offset == static_cast<int>(solver_cliques[i].size()))
                        continue;

                    // only one vertex is not removed
                    if (static_cast<int>(solver_cliques[i].size()) - 1 == offset) {
                        if (pick)
                            pick_vertex(unremoved_vertex);
                        else
                            unpick_vertex(unremoved_vertex);
                    } else
                        create_constraint(solver_cliques[i], solver_cliques[i].size() - offset);
                }
            }
        }

        if (status == sat_status::null) {
            if (flags.verbose)
                std::cerr << "! oracle fails" << std::endl;
            return false;
        }
    }

    return true;
};

void Solver::create_totalizer_from_constraint(int constraint_id) {
    assert(!constraint_literal[constraint_id]);
    assert(constraint_enabled[constraint_id]);
    assert(!is_constraint_removed(constraint_id));

    auto &constraint = constraints[constraint_id];
    assert(static_cast<int>(constraint.get_left_side().size()) > constraint.get_right_side());

    // Handle the edge case that a constraint is trivally unsatisfiable.
    // This might happen if a trivally unsatisfiable constraint is passed to create_constraint.
    if (constraint.get_right_side() < 0) {
        int always_false = sink.allocate();
        sink.emit(std::array<atom, 1>{always_false});
        constraint_literal[constraint_id] = -always_false;
        return;
    }
    assert(!constraint.get_left_side().empty());

    if (static_cast<int>(constraint.get_left_side().size()) - 1 == constraint.get_right_side()) {
        stats_.prevent_create_totalizer_of_constraint_count++;

        std::vector<int> clause;
        clause.push_back(sink.allocate());
        for (auto v : constraint.get_left_side()) {
            assert(graph_.valid_vertex(v));
            assert(x[v]);
            clause.push_back(-x[v]);
        }

        sink.emit(clause);
        constraint_literal[constraint_id] = -clause.front();
    } else {
        stats_.create_totalizer_of_constraint_count++;

        std::vector<atom> t_; // Output of the totalizer
        t_.push_back(sink.allocate());

        std::vector<atom> c_;
        for (auto v : constraint.get_left_side()) {
            assert(graph_.valid_vertex(v));
            assert(x[v]);
            c_.push_back(x[v]);
        }

        if (flags.verbose)
            std::cerr << "constructing totalizer for 0 <= q < " << constraint.get_right_side()
                      << std::endl;
        totalizer(sink, c_, t_, constraint.get_right_side(), constraint.get_right_side() + 1);
        constraint_literal[constraint_id] = -t_.front();
    }
}

Solution Solver::get_upper_bound() {
    Solution s{0};
    if (!graph_.empty_graph()) {
        if (flags.ub == ub_algorithm::greedy) {
            auto start = std::chrono::high_resolution_clock::now();
            auto vc = algorithms::upper_bound_greedy(graph_);
            auto end = std::chrono::high_resolution_clock::now();
            stats_.upper_bound_time += end - start;

            project_solution(s, vc);
            return s;
        } else if (flags.ub == ub_algorithm::approx2) {
            auto start = std::chrono::high_resolution_clock::now();
            auto vc = algorithms::upper_bound_2_guarantee(graph_);
            auto end = std::chrono::high_resolution_clock::now();
            stats_.upper_bound_time += end - start;

            project_solution(s, vc);
            return s;
        } else if (flags.ub == ub_algorithm::ils) {
            auto start = std::chrono::high_resolution_clock::now();
            auto vc = algorithms::iterated_local_search(graph_);
            auto end = std::chrono::high_resolution_clock::now();
            stats_.upper_bound_time += end - start;

            project_solution(s, vc);
            return s;
        }
    }

    std::vector<vertex_type> vc;
    for (vertex_type v : graph_.get_valid_vertices())
        vc.push_back(v);
    project_solution(s, vc);
    return s;
}

int Solver::get_lower_bound() {
    int lower_bound = 0;

    auto start = std::chrono::high_resolution_clock::now();
    if (flags.clique) {
        std::vector<std::vector<int>> cliques_;
        lower_bound = std::max(lower_bound, algorithms::clique_cover(graph_, cliques_));
    } else if (!flags.use_lp_reduction && flags.cycle)
        lower_bound = std::max(lower_bound, algorithms::cycleCoverILPRelaxation(graph_));
    auto end = std::chrono::high_resolution_clock::now();
    stats_.lower_bound_time += end - start;

    return lower_bound;
}

bool Solver::solve() {
    bool found_first_solution = false;
    int lower_bound = 0;

    decisions_.emplace_back();
    auto &rd = decisions_.back();
    rd.vertex = -1;
    rd.trail_ptr = 0;
    rd.bound = 0;

    if (encode_immediately)
        changed_vertices = INT_MAX / 2;

    // Get an initial upper bound.
    Solution ub_solution = get_upper_bound();
    stats_.initial_cover = ub_solution.size();
    best_solution = std::move(ub_solution);
    std::cerr << "initial solution: " << best_solution.size() << std::endl;

    if (clean_up()) {
        std::cerr << "instance solved by initial reductions" << std::endl;
        return true;
    }

    // Compute initial lower bound, this is necessary to invoke SAT.
    lower_bound = get_lower_bound();
    rd.bound = std::max(rd.bound, current_cardinality + lower_bound);
    std::cerr << "initial lower bound: " << rd.bound << std::endl;

    // Exclude the case that the instance is ready solved.
    if (rd.bound >= best_solution.size())
        return true;

    if (flags.only_sat) {
        search_oracle(oracle_strategy::exhaustive);

        if (flags.verbose) {
            std::cerr << "Solution improved to " << best_solution.size() << " vertices: ";
            for (size_t i = 0; i < best_solution.picked_vertices.size(); i++) {
                if (best_solution.picked_vertices[i])
                    std::cerr << i << " ";
            }
            std::cerr << std::endl;
        }

        return true;
    }

    if (flags.use_sat_solver) {
        if (search_oracle(oracle_strategy::initial_ub)) {
            std::cerr << "instance solved by initial SAT (upper bound)" << std::endl;
            return true;
        }
        if (rd.bound >= best_solution.size())
            return true;

        if (search_oracle(oracle_strategy::none)) {
            std::cerr << "instance solved by initial SAT (lower bound)" << std::endl;
            return true;
        }
        if (rd.bound >= best_solution.size())
            return true;
    }

    std::cerr << "bounds after initial SAT: [" << rd.bound << ", " << best_solution.size() << ")"
              << std::endl;

    struct sigaction new_signal, old_signal;
    new_signal.sa_handler = catch_sigxcpu;
    sigemptyset(&new_signal.sa_mask);
    new_signal.sa_flags = SA_RESTART;
    // cpu timeout
    sigaction(SIGXCPU, &new_signal, &old_signal);
    // ctrl + c
    sigaction(SIGINT, &new_signal, &old_signal);

    const double ewa_factor = 0.05;

    double reductions_ewa = graph_.num_of_vertices();
    double reductions_success_ewa = 1.0;
    double sat_discount = 1.0;

    while (!timeout_signal) {
        iteration_count++;

        if (flags.verbose) {
            std::cerr << "----------------------------" << std::endl
                      << "New iteration: " << iteration_count << std::endl
                      << "Current found solution size: " << current_cardinality << std::endl;
        }

        // Invalidate the solver if the graph changes considerably.
        if (solver_constructed_at != -1) {
            if (graph_.num_of_vertices()
                < vertices_at_construction - flags.sat_invalidation_count) {
                solver_constructed_at = -1;
                ++stats_.solver_invalidations;
            }
        }

        auto size_before_reductions = graph_.num_of_vertices();
        if (clean_up()) {
            if (!backtrack())
                return true;
            continue;
        }
        auto reductions_delta = size_before_reductions - graph_.num_of_vertices();
        bool reductions_were_effective = false;
        if (reductions_delta > 0)
            reductions_were_effective = true;

        reductions_ewa = ewa_factor * reductions_delta + (1 - ewa_factor) * reductions_ewa;
        reductions_success_ewa = ewa_factor * (reductions_delta > 10 ? 1 : 0)
                                 + (1 - ewa_factor) * reductions_success_ewa;

        if (timeout_signal)
            break;

        lower_bound = get_lower_bound();

        if (flags.use_lp_reduction) {
            int lower_bound_ = lp_reduction(lower_bound);
            if (lower_bound_ > lower_bound)
                lower_bound = lower_bound_;
        }

        auto &d = decisions_.back();
        d.bound = std::max(d.bound, current_cardinality + lower_bound);

        if (d.bound >= best_solution.size()
            || (!graph_.empty_graph() && get_next_vertex() == INT_MAX)) {
            if (!found_first_solution)
                stats_.upper_bound_count++;

            ++stats_.lower_bound_count;

            if (!backtrack())
                return true;
            continue;
        }

        if (flags.use_cc && !graph_.empty_graph()) {
            ccs_.clear();
            auto start = std::chrono::high_resolution_clock::now();
            algorithms::get_connected_components(graph_, ccs_);
            auto end = std::chrono::high_resolution_clock::now();
            stats_.cc_time += end - start;

            assert(ccs_.size() > 0);
            if (ccs_.size() > 1) {
                if (!solve_ccs())
                    return false;
                continue;
            }
        }

        // If picking all remaining vertices improves the solution, do it.
        // (This guarantees that the upper bound in the SAT part makes sense,
        // i.e., it is lower than the number of vertices in the graph).
        if (current_cardinality + graph_.num_of_vertices() < best_solution.size()) {
            if (!found_first_solution)
                found_first_solution = true;

            std::vector<vertex_type> vertices;
            for (auto v : graph_.get_valid_vertices())
                vertices.push_back(v);
            project_solution(best_solution, vertices);
            assert(best_solution.size() == current_cardinality + graph_.num_of_vertices());

            if (flags.verbose)
                std::cerr << "Solution improved to " << best_solution.size() << " vertices."
                          << std::endl;

            if (graph_.empty_graph()) {
                if (!backtrack())
                    return true;
                continue;
            }
        }

        if (flags.use_sat_solver) {
            if (use_probing && !reductions_were_effective) {
                auto start = std::chrono::high_resolution_clock::now();
                bool settled = search_oracle(oracle_strategy::probe);
                auto end = std::chrono::high_resolution_clock::now();
                stats_.probe_time += end - start;

                if (settled) {
                    ++stats_.probe_settled_count;
                    if (!backtrack())
                        return true;
                    continue;
                } else {
                    ++stats_.probe_surrendered_count;
                }
            }

            if (changed_vertices > sat_discount * flags.sat_factor * reductions_ewa) {
                changed_vertices = 0;
                auto start = std::chrono::high_resolution_clock::now();
                bool settled = search_oracle(oracle_strategy::none, 1 / reductions_success_ewa);
                auto end = std::chrono::high_resolution_clock::now();
                stats_.sat_solver_time += end - start;

                if (settled) {
                    sat_discount = 0.5;
                    ++stats_.branches_settled_count;
                    if (!backtrack())
                        return true;
                    continue;
                } else {
                    sat_discount = 1.0;
                    ++stats_.branches_surrendered_count;
                }
            }
        }

        branching();
    }

    return false;
}

void Solver::print_stats() const {
    std::cout << "branching:\n";
    if (flags.mirror_branch) {
        std::cout << "    mirror:\n";
        std::cout << "        time: " << stats_.satmir_time.count() << "\n";
        std::cout << "        num_found: " << stats_.mirror_count << "\n";
        std::cout << "        branches_successful: " << stats_.mirror_branches << std::endl;
    }
    if (flags.satellite_branch) {
        std::cout << "    satellite:\n";
        std::cout << "        time: " << stats_.satmir_time.count() << "\n";
        std::cout << "        num_found: " << stats_.satellite_count << "\n";
        std::cout << "        branches_successful: " << stats_.satellite_branches << std::endl;
    }

    std::cout << "reductions:\n";

    std::cout << "    degree0:\n";
    std::cout << "        time: " << stats_.reduction_zero_time.count() << "\n";
    std::cout << "        num_found: " << stats_.reduction_zero_count << "\n";
    std::cout << "        branches_successful: " << stats_.reduction_zero_branch_count << std::endl;

    if (flags.degree1) {
        std::cout << "    degree1:\n";
        std::cout << "        time: " << stats_.reduction_one_time.count() << "\n";
        std::cout << "        num_found: " << stats_.reduction_one_count << "\n";
        std::cout << "        branches_successful: " << stats_.reduction_one_branch_count
                  << std::endl;
    }

    if (flags.deg2) {
        std::cout << "    degree2:\n";
        std::cout << "        time: " << stats_.reduction_two_time.count() << "\n";
        std::cout << "        num_found: " << stats_.reduction_two_count << "\n";
        std::cout << "        branches_successful: " << stats_.reduction_two_branch_count
                  << std::endl;
    }

    if (flags.domination) {
        std::cout << "    domination:\n";
        std::cout << "        time: " << stats_.reduction_dom_time.count() << "\n";
        std::cout << "        num_found: " << stats_.reduction_dom_count << "\n";
        std::cout << "        branches_successful: " << stats_.reduction_dom_branch_count
                  << std::endl;
    }

    std::cout << "    simplicial:\n";
    std::cout << "        time: " << stats_.reduction_simplicial_time.count() << "\n";
    std::cout << "        num_found: " << stats_.reduction_simplicial_count << "\n";
    std::cout << "        branches_successful: " << stats_.reduction_simplicial_branch_count
              << std::endl;

    std::cout << "    unconfined:\n";
    std::cout << "        time: " << stats_.reduction_unconfined_time.count() << "\n";
    std::cout << "        num_found: " << stats_.reduction_unconfined_count << "\n";
    std::cout << "        branches_successful: " << stats_.reduction_unconfined_branch_count
              << std::endl;

    std::cout << "    funnel:\n";
    std::cout << "        time: " << stats_.reduction_funnel_time.count() << "\n";
    std::cout << "        num_found: " << stats_.reduction_funnel_count << "\n";
    std::cout << "        branches_successful: " << stats_.reduction_funnel_branch_count
              << std::endl;

    if (flags.ub != ub_algorithm::none) {
        std::cout << "upperbound:\n";
        if (flags.ub == ub_algorithm::greedy)
            std::cout << "    algo: greedy\n";
        else if (flags.ub == ub_algorithm::approx2)
            std::cout << "    algo: approx2\n";
        else if (flags.ub == ub_algorithm::ils)
            std::cout << "    algo: ils\n";
        std::cout << "    initial_cover: " << stats_.initial_cover << "\n";
        std::cout << "    time: " << stats_.upper_bound_time.count() << "\n";
        std::cout << "    branches_skipped: " << stats_.upper_bound_count << std::endl;
    }

    if (flags.clique || flags.cycle) {
        std::cout << "lowerbound:\n";
        if (flags.clique)
            std::cout << "    used_algo: clique\n";
        else if (flags.cycle)
            std::cout << "    used_algo: cycle\n";
        std::cout << "    time: " << stats_.lower_bound_time.count() << "\n";
        std::cout << "    branches_skipped: " << stats_.lower_bound_count << std::endl;
    }

    if (flags.packing) {
        std::cout << "packing:\n";
        std::cout << "    time: " << stats_.packing_time.count() << "\n";
        std::cout << "    num_constraints: " << stats_.packing_constraints << "\n";
        std::cout << "    branches_closed: " << stats_.packing_branches_closed << "\n";
        std::cout << "    reduction1_count: " << stats_.packing_reduction1_count << "\n";
        std::cout << "    reduction2_count: " << stats_.packing_reduction2_count << "\n";
        std::cout << "    reduction3_count: " << stats_.packing_reduction3_count << "\n";
        std::cout << "    rhs_zero_time: " << stats_.packing_rhs_zero_time.count() << "\n";
        std::cout << "    rhs_nonzero_time: " << stats_.packing_rhs_nonzero_time.count() << "\n";
        std::cout << "    yc_fixed_count: " << stats_.packing_yc_fixed_count << "\n";
        if (flags.constraint_totalizer) {
            std::cout << "    encode_as_clause: "
                      << stats_.prevent_create_totalizer_of_constraint_count << '\n';
            std::cout << "    encode_as_totalizer: " << stats_.create_totalizer_of_constraint_count
                      << '\n';
        }
    }

    if (flags.use_cc) {
        std::cout << "components:\n";
        std::cout << "    num_decompositions: " << stats_.cc_branches << std::endl;
        std::cout << "    time: " << stats_.cc_time.count() << "\n";
    }

    if (flags.use_sat_solver) {
        if (use_probing) {
            std::cout << "sat_probe:\n";
            std::cout << "    branches_surrendered: " << stats_.probe_surrendered_count << "\n";
            std::cout << "    branches_settled: " << stats_.probe_settled_count << "\n";
            std::cout << "    time: " << stats_.probe_time.count() << "\n";
        }
        std::cout << "sat_solver:\n";
        std::cout << "    num_constructions: " << stats_.solver_constructions << "\n";
        std::cout << "    num_invalidations: " << stats_.solver_invalidations << "\n";
        std::cout << "    num_sat: " << stats_.sat_count << "\n";
        std::cout << "    num_unsat: " << stats_.unsat_count << "\n";
        std::cout << "    branches_surrendered: " << stats_.branches_surrendered_count << "\n";
        std::cout << "    branches_settled: " << stats_.branches_settled_count << "\n";
        std::cout << "    time: " << stats_.sat_solver_time.count() << "\n";
        std::cout << "    x_fixed: " << stats_.x_fixed << '\n';
        std::cout << "    t_fixed: " << stats_.t_fixed << '\n';
    }
}

} // namespace vc_bnb
