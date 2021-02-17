#pragma once

#include <decision.hpp>
#include <vector>

struct TrailState {
    TrailState() : stack_(0) {}

    void push(action act) { stack_.push_back(act); }

    action pop() {
        auto last_decision = stack_.back();
        stack_.pop_back();

        return last_decision;
    }

    action top() { return stack_.back(); }

    action action_at(int vertex) const { return stack_[vertex]; }

    int length() const { return stack_.size(); }

private:
    std::vector<action> stack_;
};
