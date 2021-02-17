#pragma once

#include <action.hpp>
#include <chrono>
#include <constraint.hpp>
#include <decision.hpp>
#include <flags.hpp>
#include <graph.hpp>
#include <sat.hpp>
#include <solution.hpp>
#include <states.hpp>
#include <stats.hpp>
#include <trailstate.hpp>

namespace vc_bnb {

struct Solver {
    Solver(Graph &graph);

    Solver() = delete;

    bool is_picked(int vertex);
    bool is_unpicked(int vertex);
    bool is_removed(int vertex) const;
    bool is_unified(int vertex) const;
    bool got_unified(int vertex) const;
    bool is_constraint_removed(int constraint_id) const;

    std::tuple<int, int, int> get_unified_vertices(int vertex) const;

    void project_solution(Solution &out, const std::vector<vertex_type> &vc);

    int get_next_vertex();

    void branching();

    void track_decision(action act);

    bool pick_vertex(int vertex);
    bool unpick_vertex(int vertex);
    bool unpick_satellites();
    bool remove_in_constraints(int vertex, bool decrease_right_side);
    bool check_constraints();
    void remove_constraint(int constraint_id);
    void create_pick_constraint(int vertex);
    void create_unpick_constraint(int vertex);
    void create_constraint(std::vector<int> left_side, int right_side);
    void disable_constraints(int vertex);
    std::vector<Constraint> get_constraints() const { return constraints; }
    std::vector<std::vector<int>> get_vertex_to_constraints() const {
        return vertex_to_constraints;
    }

    void reduction_degree_zero();
    void reduction_degree_one();
    void reduction_degree_two();
    void reduction_domination();
    void reduction_simplicial();
    void reduction_unconfined();
    void reduction_funnel();
    int lp_reduction(int lowerbound);

    bool clean_up();

    bool solve_ccs();

    enum class oracle_strategy { none, initial_ub, exhaustive, probe };

    bool search_oracle(oracle_strategy strategy, double limit_factor = 1.0);
    void create_totalizer_from_constraint(int constraint_id);

    void undo_unified_vertex(int vertex);
    void undo_step();
    bool backtrack();

    Solution get_upper_bound();
    int get_lower_bound();

    bool solve();

    void print_solution() const;

    std::vector<bool> get_best_solution() const { return best_solution.picked_vertices; }

    void print_stats() const;

    int get_start_vertex_count() const { return start_vertex_count; }

    Flags flags;

    int solver_constructed_at = -1;
    int cardinality_at_construction;
    int lb_at_construction;
    int ub_at_construction;
    int vertices_at_construction;
    int id_limit_at_construction;

    std::vector<std::vector<int>> satellites_;
    Solution best_solution;
    int iteration_count = 0;

private:
    Graph graph_;
    TrailState trail_;
    std::vector<decision> decisions_;
    int start_vertex_count = 0;
    int changed_vertices = 0;

    std::vector<std::tuple<int, int, int>> unified_vertices_;
    states states_;

    std::vector<std::vector<int>> ccs_;
    std::vector<int> found_mirrors;
    std::vector<int> found_satellites;
    std::vector<vertex_type> alternatives;
    std::vector<std::vector<int>> solver_cliques;
    std::vector<Constraint> constraints;
    std::vector<std::vector<int>> vertex_to_constraints;
    std::vector<bool> got_unified_;
    std::vector<bool> constraint_enabled;
    std::vector<bool> constraint_removed;
    std::vector<atom> constraint_literal;
    int current_cardinality = 0;

    cadical_oracle sink;
    std::vector<atom> x; // x_v variables
    std::vector<atom> c; // Array of all x_v variables
    std::vector<atom> y; // Array of all clique variables
    std::vector<atom> t; // Output of the totalizer
    std::vector<int> new_picked_vertices;
    std::vector<int> new_unpicked_vertices;
    std::vector<int> new_constraints;
    std::vector<atom> new_folds;

    stats stats_;
};

} // namespace vc_bnb
