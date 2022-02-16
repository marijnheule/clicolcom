#ifndef __VC_REDUCTION_HH
#define __VC_REDUCTION_HH

#include <assert.h>
#include <vector>

// #include "graph.hpp"
#include "prop.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

namespace gc
{

enum class vertex_status : uint8_t {
    in_graph,
    low_degree_removed,
    indset_removed,
    dominated_removed,
    dsatur_removed,
};

template <class adjacency_struct> struct graph_reduction {
    const gc::graph<adjacency_struct>& g;
    const gc::statistics& statistics;
    std::vector<int> removed_vertices;
    std::vector<int> dominator;
    std::vector<gc::indset_constraint> constraints;
    std::vector<vertex_status> status;
    gc::bitset nodeset;
    // gc::bitset util_set;
    std::vector<int>& solution;

    explicit graph_reduction(const gc::graph<adjacency_struct>& g,
        const gc::statistics& statistics, std::vector<int>& solution)
        : g(g)
        , statistics(statistics)
        , status(g.capacity(), vertex_status::in_graph)
        , nodeset(0, g.capacity(), gc::bitset::empt)
        // , util_set(0, g.capacity(), gc::bitset::empt)
        , solution(solution)
    {
    }

    int extend_solution(const int* col, const int ub)
    {
        for (int v = 0; v < g.size(); ++v) {
            solution[g.nodes[v]] = col[v];
        }

        return extend_solution(solution, ub, true);
    }

    int extend_solution(std::vector<int>& col, const int ub, const bool full = false)
    {
				gc::bitset util_set;
        util_set.initialise(0, ub, gc::bitset::empt);

        int maxc{0};
        nodeset.copy(g.nodeset);

        for (auto v : g.nodes) {
            maxc = std::max(maxc, col[v]);
        }

        auto d{dominator.rbegin()};

        //
        // std::cout << "\n |removed| = " << removed_vertices.size() << " "
        //           << nodeset.size() << " " << nodeset << std::endl;

        for (auto i = removed_vertices.rbegin(), iend = removed_vertices.rend();
             i != iend; ++i) {

            auto v = *i;

            if (!full and status[v] != vertex_status::indset_removed)
                break;

            if (status[v] == vertex_status::dominated_removed) {
                assert(nodeset.fast_contain(*d));
                col[v] = col[*d];
                ++d;
            }

            assert(status[v] != vertex_status::in_graph);
            util_set.fill();
            for (auto u : g.matrix[v]) {
                if (!nodeset.fast_contain(u))
                    continue;
                util_set.fast_remove(col[u]);
            }

            int q{util_set.min()};

            // if (maxc < q) {
            //     std::cout << "need new color for " << v << " ("
            //               << g.matrix[v].size() << ") " << g.matrix[v]
            //               << std::endl;
            // }

            maxc = std::max(maxc, q);
            col[v] = q;

            nodeset.fast_add(v);
        }

        return maxc + 1;
    }

    template <class viterator>
    int greedy_solution(
        std::vector<int>& col, viterator first, viterator last, const int ub)
    {
				gc::bitset util_set;
      	util_set.initialise(0, ub, gc::bitset::empt);

        int maxc{0};
        nodeset.clear();

        for (auto i = first; i != last; ++i) {
            auto v = *i;

            util_set.fill();
            for (auto u : g.matrix[v]) {
                if (!nodeset.fast_contain(u))
                    continue;
                util_set.fast_remove(col[u]);
            }

            int q{util_set.min()};

            maxc = std::max(maxc, q);
            col[v] = q;

            assert(q <= ub);
            nodeset.fast_add(v);
        }

        return maxc + 1;
    }
};
}

#endif
