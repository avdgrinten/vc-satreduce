#pragma once

struct stats {
    // branching
    unsigned int mirror_count = 0;
    unsigned int mirror_branches = 0;
    unsigned int satellite_count = 0;
    unsigned int satellite_branches = 0;
    std::chrono::duration<double> satmir_time = std::chrono::duration<double>(0.0);

    // reductions
    unsigned int reduction_zero_count = 0;
    unsigned int reduction_one_count = 0;
    unsigned int reduction_two_count = 0;
    unsigned int reduction_dom_count = 0;
    unsigned int reduction_simplicial_count = 0;
    unsigned int reduction_unconfined_count = 0;
    unsigned int reduction_funnel_count = 0;
    unsigned int reduction_zero_branch_count = 0;
    unsigned int reduction_one_branch_count = 0;
    unsigned int reduction_two_branch_count = 0;
    unsigned int reduction_dom_branch_count = 0;
    unsigned int reduction_simplicial_branch_count = 0;
    unsigned int reduction_unconfined_branch_count = 0;
    unsigned int reduction_funnel_branch_count = 0;

    std::chrono::duration<double> reduction_zero_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_one_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_two_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_dom_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_simplicial_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_unconfined_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> reduction_funnel_time = std::chrono::duration<double>(0.0);

    // upper/lowerbound
    int initial_cover = 0;
    unsigned int upper_bound_count = 0;
    unsigned int lower_bound_count = 0;

    std::chrono::duration<double> upper_bound_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> lower_bound_time = std::chrono::duration<double>(0.0);

    // packing
    unsigned int packing_constraints = 0;
    unsigned int packing_branches_closed = 0;
    unsigned int packing_reduction1_count = 0;
    unsigned int packing_reduction2_count = 0;
    unsigned int packing_reduction3_count = 0;
    unsigned int packing_yc_fixed_count = 0;

    std::chrono::duration<double> packing_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> packing_rhs_zero_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> packing_rhs_nonzero_time = std::chrono::duration<double>(0.0);

    // Connected components.
    unsigned int cc_branches = 0;
    std::chrono::duration<double> cc_time = std::chrono::duration<double>(0.0);

    // sat solver
    unsigned int solver_constructions = 0;
    unsigned int solver_invalidations = 0;
    unsigned int sat_count = 0;
    unsigned int unsat_count = 0;
    unsigned int probe_settled_count = 0;
    unsigned int probe_surrendered_count = 0;
    unsigned int branches_settled_count = 0;
    unsigned int branches_surrendered_count = 0;
    unsigned int create_totalizer_of_constraint_count = 0;
    unsigned int prevent_create_totalizer_of_constraint_count = 0;
    unsigned int x_fixed = 0;
    unsigned int t_fixed = 0;

    std::chrono::duration<double> probe_time = std::chrono::duration<double>(0.0);
    std::chrono::duration<double> sat_solver_time = std::chrono::duration<double>(0.0);

    stats &operator+(const stats &other) {
        this->mirror_count = +other.mirror_count;
        this->mirror_branches += other.mirror_branches;
        this->satellite_count += other.satellite_count;
        this->satellite_branches += other.satellite_branches;
        this->satmir_time += other.satmir_time;

        this->reduction_zero_count += other.reduction_zero_count;
        this->reduction_one_count += other.reduction_one_count;
        this->reduction_two_count += other.reduction_two_count;
        this->reduction_dom_count += other.reduction_dom_count;
        this->reduction_simplicial_count += other.reduction_simplicial_count;
        this->reduction_unconfined_count += other.reduction_unconfined_count;
        this->reduction_funnel_count += other.reduction_funnel_count;
        this->reduction_zero_branch_count += other.reduction_zero_branch_count;
        this->reduction_one_branch_count += other.reduction_one_branch_count;
        this->reduction_two_branch_count += other.reduction_two_branch_count;
        this->reduction_dom_branch_count += other.reduction_dom_branch_count;
        this->reduction_simplicial_branch_count += other.reduction_simplicial_branch_count;
        this->reduction_unconfined_branch_count += other.reduction_unconfined_branch_count;
        this->reduction_funnel_branch_count += other.reduction_funnel_branch_count;

        this->reduction_zero_time += other.reduction_zero_time;
        this->reduction_one_time += other.reduction_one_time;
        this->reduction_two_time += other.reduction_two_time;
        this->reduction_dom_time += other.reduction_dom_time;
        this->reduction_simplicial_time += other.reduction_simplicial_time;
        this->reduction_unconfined_time += other.reduction_unconfined_time;
        this->reduction_funnel_time += other.reduction_funnel_time;

        this->initial_cover += other.initial_cover;
        this->upper_bound_count += other.upper_bound_count;
        this->lower_bound_count += other.lower_bound_count;

        this->upper_bound_time += other.upper_bound_time;
        this->lower_bound_time += other.lower_bound_time;

        this->packing_constraints += other.packing_constraints;
        this->packing_branches_closed += other.packing_branches_closed;
        this->packing_reduction1_count += other.packing_reduction1_count;
        this->packing_reduction2_count += other.packing_reduction2_count;
        this->packing_reduction3_count += other.packing_reduction3_count;
        this->packing_yc_fixed_count += other.packing_yc_fixed_count;

        this->packing_time += other.packing_time;
        this->packing_rhs_zero_time += other.packing_rhs_zero_time;
        this->packing_rhs_nonzero_time += other.packing_rhs_nonzero_time;

        this->cc_branches += other.cc_branches;
        this->cc_time += other.cc_time;

        this->solver_constructions += other.solver_constructions;
        this->solver_invalidations += other.solver_invalidations;
        this->sat_count += other.sat_count;
        this->unsat_count += other.unsat_count;
        this->probe_settled_count += other.probe_settled_count;
        this->probe_surrendered_count += other.probe_surrendered_count;
        this->branches_settled_count += other.branches_settled_count;
        this->branches_surrendered_count += other.branches_surrendered_count;
        this->create_totalizer_of_constraint_count += other.create_totalizer_of_constraint_count;
        this->prevent_create_totalizer_of_constraint_count +=
            other.prevent_create_totalizer_of_constraint_count;
        this->x_fixed += other.x_fixed;
        this->t_fixed += other.t_fixed;

        this->probe_time += other.probe_time;
        this->sat_solver_time += other.sat_solver_time;

        return *this;
    }
};
