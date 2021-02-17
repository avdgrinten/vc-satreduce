#pragma once

#include <action.hpp>
#include <graph.hpp>

enum class alternative_type { none, include, exclude, exclude_satellites };

struct decision {
    vertex_type vertex;
    int trail_ptr;
    alternative_type alternative;
    int bound;
};

struct action {
    action(vc_bnb::action_type type) : type(type) { constraint_id = -1; }

    action(int vertex, vc_bnb::action_type type) : vertex(vertex), type(type) {
        constraint_id = -1;
    };

    action(int vertex, vc_bnb::action_type type, int constraint_id)
        : vertex(vertex), type(type), constraint_id(constraint_id){};

    action() : vertex(-1), type(vc_bnb::action_type::invalid) {}

    void set_action(vc_bnb::action_type action_) { type = action_; };

    int vertex;
    vc_bnb::action_type type;
    int central;
    int contracted[2];
    int constraint_id;
};
