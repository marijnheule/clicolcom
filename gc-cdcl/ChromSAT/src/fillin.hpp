#ifndef __GC_FINLLIN_HPP
#define __GC_FINLLIN_HPP

#include "graph.hpp"

namespace gc
{

template <class graph_struct> struct minfill_buffer {

    const graph_struct& orig;

    // result
    graph_struct g;
    std::vector<int> order;
    std::vector<std::pair<int, int>> fillin;
    int width{0};

    // buffers
    std::vector<int> degrees;
    std::vector<std::list<int>::iterator> iterators;
    bitset ordered;
    std::vector<std::list<int>> buckets;
    bitset util_set;

    void change_degree(int u, int diff)
    {
        auto& ud = degrees[u];
        buckets[ud].erase(iterators[u]);
        ud += diff;
        if (ud >= static_cast<int>(buckets.size()))
            buckets.resize(ud + 1);
        buckets[ud].push_front(u);
        iterators[u] = buckets[ud].begin();
    }

    explicit minfill_buffer(const graph_struct& g)
        : orig(g)
        , g(g)
        , degrees(g.capacity())
        , iterators(g.capacity())
        , ordered(0, g.capacity() - 1, bitset::empt)
        , util_set(0, g.capacity() - 1, bitset::empt)
    {
    }

    void minfill()
    {
        for (auto v : g.nodes) {
            util_set.copy(g.matrix[v]);
            util_set.intersect_with(g.nodeset);
            auto vd = util_set.size();
            if (vd >= buckets.size())
                buckets.resize(vd + 1);
            buckets[vd].push_front(v);
            degrees[v] = vd;
            iterators[v] = buckets[vd].begin();
            ordered.clear();
        }

        while (true) {
            size_t i{0};
            for (; i != buckets.size(); ++i)
                if (!buckets[i].empty())
                    break;
            if (i == buckets.size())
                break;
            // Here i is the first non-empty buckets = the lowest degree
            width = std::max(width, static_cast<int>(i));
            // So is width
            auto v = buckets[i].back();
            order.push_back(v);
            buckets[i].pop_back();
            ordered.fast_add(v);
            for (auto u : g.matrix[v]) {
                if (ordered.fast_contain(u) || !g.nodeset.fast_contain(u))
                    continue;
                util_set.copy(g.matrix[v]);
                util_set.setminus_with(g.matrix[u]);
                util_set.setminus_with(ordered);
                util_set.intersect_with(g.nodeset);
                for (auto w : util_set) {
                    if (w == u)
                        continue;
                    assert(!g.matrix[w].fast_contain(u));
                    g.add_edge(u, w);
                    fillin.push_back(std::pair<int, int>{u, w});
                    change_degree(u, 1);
                    change_degree(w, 1);
                }

                change_degree(u, -1);
            }
        }
    }
};

} // namespace gc

#endif
