#pragma once

#include <graph.hpp>
#include <vector>

namespace algorithms {

// Returns (number of vertices - number of cliques).
int clique_cover(const Graph &g, std::vector<std::vector<int>> &cliques);

void get_isolated_vertices(Graph &g, std::vector<int> &isolted_vertives);

void find_satellites(const Graph &g, int v, std::vector<int> &satellites);

void find_mirrors(const Graph &g, int v, std::vector<int> &mirrors);

void find_mirrors_and_satellites(const Graph &graph, int v, std::vector<int> &mirrors,
                                 std::vector<int> &satellites);

std::vector<int> upper_bound_2_guarantee(const Graph &g);

std::vector<int> upper_bound_greedy(const Graph &g);

std::vector<int> iterated_local_search(const Graph &g);

void get_connected_components(const Graph &g, std::vector<std::vector<int>> &cc);

void hopcroftKarp(const Graph &G, std::vector<int> &, std::vector<int> &);

void buildCover(std::vector<int> &, std::vector<int> &, const Graph &G, const std::vector<int> &,
                const std::vector<int> &);

int cycleCoverILPRelaxation(const Graph &G);
int cycleCoverILPRelaxation(const Graph &G, std::vector<int> &add, std::vector<int> &remove);

} // namespace algorithms
