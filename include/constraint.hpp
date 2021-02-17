#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include <action.hpp>

struct Constraint {
    Constraint(std::vector<int> lhs, int rhs)
        : original_vertices(std::move(lhs)), original_right_side(rhs) {
        vertices = original_vertices;
        right_side = original_right_side;
    }

    void update(int vertex, bool decrease_right_side) {
        bool found = false;
        for (size_t i = 0; i < vertices.size(); i++) {
            auto xv = vertices[i];

            if (xv == vertex) {
                found = true;
                vertices.erase(vertices.begin() + i);
                break;
            }
        }
        assert(found);

        if (decrease_right_side)
            right_side--;
    }

    void backtrack(int vertex, bool increase_right_side) {
        for (auto v : vertices)
            assert(v != vertex);

        vertices.push_back(vertex);
        std::sort(vertices.begin(), vertices.end());

        if (increase_right_side)
            right_side++;
    }

    const std::vector<int> &get_original_left_side() const { return original_vertices; }

    const std::vector<int> &get_left_side() const { return vertices; }

    int get_right_side() const { return right_side; }

    bool equals(const std::vector<int> &left_side, int right_side) const {
        return this->right_side == right_side && this->vertices == left_side;
    }

private:
    std::vector<int> original_vertices;
    std::vector<int> vertices;
    int original_right_side;
    int right_side;
};
