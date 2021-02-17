#pragma once

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include <cadical.hpp>

using atom = int;

enum class sat_status { null, sat, unsat };

struct cnf_file {
    cnf_file(std::string path) : path_{std::move(path)} {}

    cnf_file(const cnf_file &) = delete;

    cnf_file &operator=(const cnf_file &) = delete;

    atom allocate() { return ++atom_counter_; }

    template <typename Range>
    void emit(Range c) {
        for (atom x : c)
            encoding_ << x << ' ';
        encoding_ << "0\n";
        num_clauses_++;
    }

    void comment(const char *s) { encoding_ << "c " << s << '\n'; }

    void finalize() {
        // TODO: check for I/O errors.
        std::ofstream out{path_};
        out << "p cnf " << atom_counter_ << " " << num_clauses_ << "\n";
        out << encoding_.str();
    }

private:
    std::string path_;
    std::stringstream encoding_;
    int atom_counter_ = 0;
    uint64_t num_clauses_ = 0;
};

struct cadical_oracle {
    cadical_oracle() { reset(); }

    atom allocate() { return ++atom_counter_; }

    void reset() {
        atom_counter_ = 0;
        solver_ = std::unique_ptr<CaDiCaL::Solver>(new CaDiCaL::Solver);
    }

    unsigned int num_atoms() { return atom_counter_; }

    template <typename Range>
    void emit(Range c) {
        assert(c.size() > 0); // prevents 'empty original clause'
        for (atom x : c) {
            assert(x != 0); // prevents 'empty original clause'
            solver_->add(x);
        }
        solver_->add(0);
    }

    void assume(atom x) { solver_->assume(x); }

    int value(atom x) { return solver_->val(x); }

    int fixed(atom x) { return solver_->fixed(x); }

    void configure_sat() {
        auto success = solver_->configure("sat");
        assert(success && "configuration does not exist");
    }

    void configure_unsat() {
        auto success = solver_->configure("unsat");
        assert(success && "configuration does not exist");
    }

    sat_status solve(int limit) {
        solver_->limit("conflicts", limit);
        int status = solver_->solve();
        if (print_stats_)
            solver_->statistics();
        if (status == 10) {
            return sat_status::sat;
        } else if (status == 20) {
            return sat_status::unsat;
        } else if (!status) {
            return sat_status::null;
        } else {
            std::cerr << "unexpected result from SAT solver: " << status << std::endl;
            abort();
        }
    }

    void enable_stats(bool p) { print_stats_ = p; }

private:
    std::unique_ptr<CaDiCaL::Solver> solver_;
    int atom_counter_;
    bool print_stats_ = false;
};

// Construct a totalizer that can check whether: sum_{i \in in} x_i <= q for any lb <= q < ub.
// The totalizer assumes (but does not enforce):
//     sum x_i >= lb
//     sum x_i < ub
// The "output" t satisfies:
//     t[0] <=> sum x_i >= lb + 1
//     t[1] <=> sum x_i >= lb + 2
//     ...
//     t[ub - lb - 2] <=> sum x_i >= lb + 1 + (ub - lb - 2) = ub - 1
//     t[ub - lb - 1] <=> sum x_i >= lb + 1 + (ub - lb - 1) = ub
template <typename S>
void totalizer(S &sink, std::vector<atom> in, std::vector<atom> t, int lb, int ub) {
    assert(static_cast<int>(in.size()) >= 1);
    assert(lb < ub);
    assert(static_cast<int>(t.size()) == ub - lb);
    assert(lb >= 0);
    assert(ub <= static_cast<int>(in.size()));

    if (in.size() == 1) {
        assert(lb == 0);
        assert(ub == 1);

        std::vector<atom> c;
        c.push_back(-in[0]);
        c.push_back(t[0]);
        sink.emit(c);
        return;
    }

    std::vector<std::tuple<std::vector<atom>, std::vector<atom>, int, int>> stack;
    stack.push_back({std::move(in), t, lb, ub});
    while (!stack.empty()) {
        std::tie(in, t, lb, ub) = std::move(stack.back());
        stack.pop_back();
        assert(in.size() >= 2);

        if (in.size() == 2) {
            if (lb == 0 && ub == 1) {
                std::vector<atom> c;
                c.push_back(-in[0]);
                c.push_back(t[0]);
                sink.emit(c);

                c.clear();
                c.push_back(-in[1]);
                c.push_back(t[0]);
                sink.emit(c);
            } else {
                assert(lb == 0); // TODO: Handle the lb == 1 && ub == 2 case.
                assert(ub == 2);

                std::vector<atom> c;
                c.push_back(-in[0]);
                c.push_back(-in[1]);
                c.push_back(t[1]);
                sink.emit(c);

                c.clear();
                c.push_back(-in[0]);
                c.push_back(t[0]);
                sink.emit(c);

                c.clear();
                c.push_back(-in[1]);
                c.push_back(t[0]);
                sink.emit(c);
            }
            continue;
        }

        // split t into children l_t and r_t
        std::vector<atom> l_t;
        std::vector<atom> r_t;

        std::vector<atom> l_in{in.begin(), in.begin() + in.size() / 2};
        std::vector<atom> r_in{in.begin() + in.size() / 2, in.end()};
        assert(l_in.size() >= 1);
        assert(r_in.size() >= 1);

        // TODO: choose the LB and UB in a smart way to make the totalizer smaller.
        int l_lb = std::max(lb - static_cast<int>(r_in.size()), 0);
        int l_ub = std::min(l_lb + static_cast<int>(l_in.size()), ub);
        int r_lb = std::max(lb - static_cast<int>(l_in.size()), 0);
        int r_ub = std::min(r_lb + static_cast<int>(r_in.size()), ub);
        assert(l_lb < l_ub);
        assert(r_lb < r_ub);

        if (l_in.size() == 1) {
            assert(l_lb == 0);
            assert(l_ub == 1);
            l_t.push_back(l_in[0]);
        } else {
            for (int i = 0; i < l_ub - l_lb; i++)
                l_t.push_back(sink.allocate()); // allocate extra t
            stack.push_back({std::move(l_in), l_t, l_lb, l_ub});
        }

        if (r_in.size() == 1) {
            assert(r_lb == 0);
            assert(r_ub == 1);
            r_t.push_back(r_in[0]);
        } else {
            for (int j = 0; j < r_ub - r_lb; j++)
                r_t.push_back(sink.allocate()); // allocate extra t
            stack.push_back({std::move(r_in), r_t, r_lb, r_ub});
        }

        // Note that cases (*) below imply that the lower bound is not satisfied.
        // We could add clauses to force the input to be 1 in these cases,
        // but they are generally not useful (as they will not appear in SAT conflicts).

        for (size_t i = 0; i < l_t.size(); i++) {
            int l_v = l_lb + i + 1;
            if (l_v < lb + 1) // (*)
                continue;
            assert(l_v <= ub);
            std::vector<atom> c;
            c.push_back(-l_t[i]);
            c.push_back(t[l_v - lb - 1]);
            sink.emit(c);
        }

        for (size_t j = 0; j < r_t.size(); j++) {
            int r_v = r_lb + j + 1;
            if (r_v < lb + 1) // (*)
                continue;
            assert(r_v <= ub);
            std::vector<atom> c;
            c.push_back(-r_t[j]);
            c.push_back(t[r_v - lb - 1]);
            sink.emit(c);
        }

        for (size_t i = 0; i < l_t.size(); i++) {
            for (size_t j = 0; j < r_t.size(); j++) {
                int l_v = l_lb + i + 1;
                int r_v = r_lb + j + 1;
                if (l_v + r_v < lb + 1) // (*)
                    continue;
                if (l_v + r_v > ub)
                    continue;
                std::vector<atom> c;
                c.push_back(-l_t[i]);
                c.push_back(-r_t[j]);
                c.push_back(t[l_v + r_v - lb - 1]);
                sink.emit(c);
            }
        }
    }
}
