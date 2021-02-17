#pragma once

#include <action.hpp>
#include <decision.hpp>

struct states {
    states(int size) { array_.resize(size, vc_bnb::action_type::invalid); }

    vc_bnb::action_type get_state(int vertex) {
        assert(vertex < static_cast<int>(array_.size()));
        return array_[vertex];
    }

    void update_state(int vertex, vc_bnb::action_type type) {
        assert(type != vc_bnb::action_type::invalid);
        assert(vertex < static_cast<int>(array_.size()));
        assert(array_[vertex] == vc_bnb::action_type::invalid);
        array_[vertex] = type;
    }

    void reset_state(int vertex) {
        assert(vertex < static_cast<int>(array_.size()));
        assert(array_[vertex] != vc_bnb::action_type::invalid);
        array_[vertex] = vc_bnb::action_type::invalid;
    }

    void add(int vertex) {
        assert(vertex == static_cast<int>(array_.size()));
        array_.push_back(vc_bnb::action_type::invalid);
    }

    void remove(int vertex) {
        assert(vertex == static_cast<int>(array_.size()) - 1);
        assert(array_[vertex] == vc_bnb::action_type::invalid);
        array_.pop_back();
    }

private:
    std::vector<vc_bnb::action_type> array_;
};
