#ifndef __GC_VC_HPP
#define __GC_VC_HPP

#include "graph.hpp"

namespace gc {

// TODO: we could use an abstract base class with a find_is() method
// to try various methods of finding good ISes, for example
// vertex-cover-cdcl or local search

template <typename graph_type>
class degeneracy_vc_solver
{
public:
    const graph_type& g;
    std::vector<int> order;
    std::vector<int> degrees;
    std::vector<std::list<int>::iterator> iterators;
    std::vector<bool> ordered;
    std::vector<std::list<int>> buckets;

    int core;
public:
    degeneracy_vc_solver(const graph_type& g)
        : g(g)
        , degrees(g.capacity())
        , iterators(g.capacity())
        , ordered(g.capacity())
    {
    }

    // find degeneracy on the complement of g. return the position of
    // the first vertex in the order so that all remaining vertices
    // have 0 degree
    int compl_degeneracy();

    // return a bitset containing an independent set of g
    bitset find_is();
};

template <typename graph_type>
int degeneracy_vc_solver<graph_type>::compl_degeneracy()
{
    for (auto v : this->g.nodes) {
        auto vd = this->g.matrix[v].size();
        if (vd >= buckets.size())
            buckets.resize(vd + 1);
        buckets[vd].push_front(v);
        degrees[v] = vd;
        iterators[v] = buckets[vd].begin();
        ordered[v] = false;
    }
		
    while (true) {
        int i = buckets.size()-1;
        for (; i >= 0; --i)
            if (!buckets[i].empty())
                break;
        if (i < 0)
            break;

        if (i == 0) {
            int rv = order.size();
            for (auto u : buckets[0])
                order.push_back(u);
            return rv;
        }

        auto v = buckets[i].back();
        order.push_back(v);
        buckets[i].pop_back();
        ordered[v] = true;
        for (auto u : this->g.matrix[v]) {
            if (ordered[u])
                continue;
            auto& ud = degrees[u];
            buckets[ud].erase(iterators[u]);
            --ud;
            buckets[ud].push_front(u);
            iterators[u] = buckets[ud].begin();
        }
    }
    // unreachable
    assert(0);
    return order.size()-1;
}

template <typename graph_type>
bitset degeneracy_vc_solver<graph_type>::find_is()
{
    int idx = compl_degeneracy();
    bitset bs(0, g.capacity(), bitset::empt);

    for (int i = idx, iend = static_cast<int>(order.size()); i != iend; ++i)
        bs.fast_add(order[i]);
    return bs;
}

}

#endif
