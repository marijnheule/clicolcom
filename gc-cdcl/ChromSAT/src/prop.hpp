#ifndef __GC_PROP_HPP
#define __GC_PROP_HPP

#include "graph.hpp"
#include "minicsp/core/solver.hpp"
#include "options.hpp"
#include "varmap.hpp"

#include <boost/optional.hpp>

namespace gc
{

using dense_graph = graph<bitset>;

struct statistics;

// constraint on a subset of the vertices of the graph that get
// added during preprocessing. We place a tighter upper bound on
// these sets of vertices
struct indset_constraint {

    template <class adjacency_struct>
    indset_constraint(adjacency_struct& scope, int s)
        // : bvs(0, scope.max(), bitset::empt)
        : source(s)
    {
        vs.reserve(scope.size());
        for (auto v : scope)
            vs.push_back(v);
        std::sort(begin(vs), end(vs));
    }

    // the vertices in the local constraint
    std::vector<int> vs;
    // where it came from (i.e., which vertex's neighborhood is
    // it). for debugging
    int source;
};

struct fillin_info {
    std::vector<std::pair<int, int>> edges;
    // vertices, in order
    std::vector<int> order;
    // reverse lookup for order: order of each vertex
    std::vector<int> revorder;

    fillin_info() = default;
    fillin_info(
        std::vector<std::pair<int, int>>&& edges, std::vector<int>&& order)
        : edges(edges)
        , order(order)
    {
        revorder.resize(order.size());
        for (size_t i = 0, iend = order.size(); i != iend; ++i)
            revorder[order[i]] = i;
    }
    fillin_info(const fillin_info&) = default;
    fillin_info(fillin_info&&) = default;
    fillin_info& operator=(const fillin_info&) = default;
    fillin_info& operator=(fillin_info&&) = default;
};

struct cons_base {
    minicsp::Solver& s;
    const options& opt;
    dense_graph& g;
    dense_graph fg;

    int ub; // ub of the original graph. [i.e., counting the vertices of the IS]
    int bestlb{0};
    clique_finder<bitset> cf; // for cliques
    clique_finder<bitset> ccf; // for clique covers
    minicsp::backtrackable<int> lastlb;

    struct varinfo_t {
        int u{-1}, v{-1};

        inline bool empty() const { return u < 0; }
    };
    // indexed by varid
    std::vector<varinfo_t> varinfo;

    // for incremental clique finding. must be here because it is
    // reset when backtracking
    bool saved_cliques{false};

    dense_graph create_filled_graph(boost::optional<fillin_info> fillin);

    explicit cons_base(minicsp::Solver& s, const options& opt, dense_graph& g,
        boost::optional<fillin_info> fillin)
        : s(s)
        , opt(opt)
        , g(g)
        , fg(create_filled_graph(fillin))
        , cf(g, (opt.memlimit < 0 ? g.capacity()
                                  : std::max(10, 100000 / g.capacity())))
        , ccf(g, (opt.indset_lb ? g.capacity() : 10))
        , lastlb(s)
        , lastdlvl(s)
    {
    }

    void sync_graph()
    {
        if (*lastdlvl < g.current_checkpoint()) {
            g.restore(*lastdlvl);
            *lastdlvl = g.current_checkpoint();
            saved_cliques = false;
        }
        while (s.decisionLevel() > g.current_checkpoint()) {
            g.checkpoint();
            *lastdlvl = g.current_checkpoint();
        }
#ifdef UPDATE_FG
        if (opt.fillin) {
            if (*lastdlvl < fg.current_checkpoint()) {
                fg.restore(*lastdlvl);
                assert(*lastdlvl == fg.current_checkpoint());
            }
            while (s.decisionLevel() > fg.current_checkpoint()) {
                fg.checkpoint();
                assert(*lastdlvl = fg.current_checkpoint());
            }
        }
#endif
    }

protected:
    // something which is kept in sync by the solver. if it diverges
    // from what the graph thinks, it means we have
    // backtracked/restarted, so we should resync to that point
    minicsp::backtrackable<int> lastdlvl;
};

cons_base* post_gc_constraint(minicsp::Solver& s, dense_graph& g,
    boost::optional<fillin_info> fillin, const varmap& vars,
    const std::vector<indset_constraint>& isconses,
    const std::vector<int>& vertex_map, const options& opt, statistics& stat);

} // namespace gc

#endif
