#pragma once

#include <boost/program_options.hpp>

namespace po = boost::program_options;

enum class ub_algorithm { none, greedy, approx2, ils };

struct Flags {
    bool verbose = false;
    bool satellite_branch = true;
    bool mirror_branch = true;
    ub_algorithm ub = ub_algorithm::ils;
    bool clique = true;
    bool deg2 = true;
    bool use_lp_reduction = false;
    bool use_cc = true;
    bool domination = true;
    bool cycle = true;
    bool degree1 = true;
    bool packing = false;
    bool use_funnel = false;
    bool unpick_constraint = true;
    bool pick_constraint = true;
    bool use_sat_solver = true;
    bool vertex_encoding = false;
    bool constraint_totalizer = true;
    int sat_conflict_limit = 8192;
    int sat_invalidation_count = 64;
    double sat_factor = 3;
    bool only_sat = false;
    bool with_reductions = true;
    bool ub_to_lb = false;
    bool check_fixed_yc = false;

    Flags() = default;

    Flags(po::variables_map vm) {
        if (vm.count("verbose"))
            verbose = true;

        if (vm.count("no-initial-solution"))
            ub = ub_algorithm::none;

        if (vm.count("lp-reduction"))
            use_lp_reduction = true;
        if (vm.count("packing"))
            packing = true;
        if (vm.count("funnel"))
            use_funnel = true;

        if (vm.count("no-sat"))
            use_sat_solver = false;

        if (vm.count("vertex-encoding"))
            vertex_encoding = true;

        if (vm.count("constraint_totalizer"))
            constraint_totalizer = true;

        if (vm.count("sat-conflict-limit"))
            sat_conflict_limit = vm["sat-conflict-limit"].as<int>();
        if (vm.count("sat-invalidation-count"))
            sat_invalidation_count = vm["sat-invalidation-count"].as<int>();
        if (vm.count("sat-factor"))
            sat_factor = vm["sat-factor"].as<double>();

        if (vm.count("only-sat"))
            only_sat = true;

        if (vm.count("with-reductions"))
            with_reductions = true;

        if (vm.count("ub-to-lb"))
            ub_to_lb = true;

        if (vm.count("check-fixed-yc"))
            check_fixed_yc = true;
    }
};
