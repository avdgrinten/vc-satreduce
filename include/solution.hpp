#pragma once

#include <vector>

struct Solution {
    Solution(int vertex_count) : picked_vertices(vertex_count, false) {}

    std::vector<bool> picked_vertices;

    void select_vertex(int vertex) {
        assert(vertex < static_cast<int>(picked_vertices.size()));
        assert(picked_vertices[vertex] == false);

        picked_vertices[vertex] = true;
    }

    void unselect_vertex(int vertex) {
        assert(vertex < static_cast<int>(picked_vertices.size()));
        assert(picked_vertices[vertex] == true);

        picked_vertices[vertex] = false;
    }

    int size() const {
        int count = 0;

        for (size_t i = 0; i < picked_vertices.size(); i++)
            if (picked_vertices[i])
                count++;

        return count;
    }
};
