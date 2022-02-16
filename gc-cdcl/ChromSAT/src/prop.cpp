#include "prop.hpp"
#include "minicsp/core/utils.hpp"
#include "utils.hpp"
#include "statistics.hpp"
#include "mycielski.hpp"

#include <cmath>

#include <iomanip>



// #define NICE_TRACE
// #define DEBUG_IS
// #define DEBUG_CC
#ifdef DEBUG_CC
#define CCOUT                                                                  \
    if (0)                                                                     \
        ;                                                                      \
    else                                                                       \
        std::cout
#else
#define CCOUT                                                                  \
    if (1)                                                                     \
        ;                                                                      \
    else                                                                       \
        std::cout
#endif

namespace gc
{

using namespace minicsp;

dense_graph cons_base::create_filled_graph(boost::optional<fillin_info> fillin)
{
    assert((fillin && opt.fillin) || (!fillin && !opt.fillin));
    if (opt.fillin) {
        dense_graph filled{g};
        for (auto e : fillin->edges)
            filled.add_edge(e.first, e.second);
        return filled;
    } else {
        dense_graph g;
        return g;
    }
}

struct bfs_state {
    Solver& s;
    const dense_graph& g;
    const dense_graph& fg;
    bitset upart, vpart, allpart;
    bitset util_set;
    enum vertex_color { WHITE, GREY, BLACK };
    std::vector<vertex_color> col;
    std::vector<int> uparent, vparent, parent;
    int rootu, rootv;

    bool single_partition;

    std::list<int> Q;

    template <typename F>
    auto visit_edge(int w, int x, F f) -> decltype(f(w, x))
    {
        if (x != rootu && x != rootv && !single_partition) {
            assert(vparent[x] != -1);
            assert(uparent[x] != -1);
            // assert(fg.origmatrix[x].fast_contain(vparent[x]));
            // assert(fg.origmatrix[x].fast_contain(uparent[x]));
            // assert(fg.origmatrix[vparent[x]].fast_contain(uparent[x]));
        }
        // std::cout << "BFS visiting edge " << w << "->" << x << "\n";
        return f(w, x);
    }

    template <typename F> auto do_bfs(int w, F f) -> decltype(f(0, 0))
    {
        typedef decltype(f(0, 0)) rettype;
        rettype rv{};
        util_set.copy(fg.origmatrix[w]);
        util_set.intersect_with(allpart);
        bool inv = vpart.fast_contain(w);
        for (auto x : util_set) {
            if (col[x] == BLACK)
                continue;
            if (col[x] == GREY) {
                rv = visit_edge(w, x, f);
                if (rv)
                    return rv;
                continue;
            }
            Q.push_back(x);
            col[x] = GREY;
            parent[x] = w;
            if (inv) {
                vparent[x] = w;
                uparent[x] = uparent[w];
            } else {
                vparent[x] = vparent[w];
                uparent[x] = w;
            }
            rv = visit_edge(w, x, f);
            if (rv)
                return rv;
        }
        col[w] = BLACK;
        return rv;
    }

    template <typename F> auto do_bfs(F f) -> decltype(f(0, 0))
    {
        typedef decltype(f(0, 0)) rettype;
        rettype rv{};
        while(!Q.empty()) {
            int v = Q.front();
            // std::cout << "BFS visiting " << v << "\n";
            Q.pop_front();
            rv = do_bfs(v, f);
            if (rv || Q.empty())
                return rv;
        }
        return rv;
    }

    template <typename F> auto bfs_inner(int u, int v, F f) -> decltype(f(0, 0))
    {
        upart.clear();
        for (auto up : g.partition[g.rep_of[u]]) {
            col[up] = WHITE;
            upart.fast_add(up);
            uparent[up] = -1;
            vparent[up] = -1;
            parent[up] = -1;
        }
        // std::cout << "upart = " << upart << "\n";
        vpart.clear();
        for (auto vp : g.partition[g.rep_of[v]]) {
            col[vp] = WHITE;
            vpart.fast_add(vp);
            uparent[vp] = -1;
            vparent[vp] = -1;
            parent[vp] = -1;
        }
        // std::cout << "vpart = " << vpart << "\n";
        allpart.copy(upart);
        allpart.union_with(vpart);
        Q.clear();
        Q.push_back(v);
        uparent[v] = u;
        vparent[u] = v;
        rootu = u;
        rootv = v;
        return do_bfs(f);
    }

    template <typename F> auto bfs(int u, int v, F f) -> decltype(f(0, 0))
    {
        single_partition = false;
        return bfs_inner(u, v, f);
    }

    template <typename F> auto bfs(int u, F f) -> decltype(f(0, 0))
    {
        single_partition = true;
        return bfs_inner(u, u, f);
    }

    bfs_state(Solver& s, const dense_graph& g, const dense_graph& fg)
        : s(s)
        , g(g)
        , fg(fg)
        , upart(0, g.capacity() - 1, bitset::empt)
        , vpart(0, g.capacity() - 1, bitset::empt)
        , allpart(0, g.capacity() - 1, bitset::empt)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , col(g.capacity())
        , uparent(g.capacity())
        , vparent(g.capacity())
        , parent(g.capacity())
    {
    }
};

class gc_constraint : public minicsp::cons, public cons_base
{
public:
    mycielskan_subgraph_finder<bitset> mf;

    // used in wake() when using fillin
    bfs_state bfs;
    // used in merge_fillin to add edges by transition in its second
    // phase
    std::vector<std::pair<int, int>> bfsroots;
    bitset needbfs;

    // a buffer for merge_fillin
    bitset vpartcopy;

    const varmap& vars;
    std::vector<indset_constraint> isconses;
    std::vector<gc::bitset> isrep;

    // used with learning=DYNAMIC_REPS
    std::vector<int> dynamic_reps;

    statistics& stat;

    // for the adaptive bound policy: this is set to true on conflicts
    // and reset to false after the stronger policy runs
    bool run_expensive_bound{false};
    minicsp::Solver::clause_callback_t adaptive_callback;
    options::dual_policy bound_source;

    // the usual assortment of vecs, vectors and bitsets to avoid
    // reallocations
    vec<Lit> reason;
    std::vector<int> culprit;
    bitset diffuv, diffvu;

    bitset util_set, util_set2;
    bitset neighborhood;
    bitset ordering_tmp, ordering_forbidden, ordering_removed_bs;

    // the representative chosen for each partition during explanation
    std::vector<int> expl_reps;
    // for fillin: if the rep is missing some edges, we need to use
    // more than one representative per partition
    std::vector<std::vector<int>> expl_reps_extra;

    std::vector<int> ordering_removed;

    std::vector<int> heuristic;

    // for incremental clique finding
    std::vector<int> altered_vertices;
    bitset altered_vertex_set;

    // for the pruning thing
    std::vector<int> changed_vertices;
    std::vector<arc> changed_edges;
    // std::vector<std::vector<int>> cliques_of;

public:
    gc_constraint(Solver& solver, dense_graph& g,
        boost::optional<fillin_info> fillin, const varmap& tvars,
        const std::vector<indset_constraint>& isconses, const options& opt,
        statistics& stat)
        : cons_base(solver, opt, g, fillin)
        , mf(g, cf, false)
        , bfs(s, g, fg)
        , bfsroots(g.capacity())
        , needbfs(0, g.capacity() - 1, bitset::empt)
        , vpartcopy(0, g.capacity() - 1, bitset::empt)
        , vars(tvars)
        , isconses(isconses)
        , isrep(isconses.size())
        , dynamic_reps(g.capacity())
        , stat(stat)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , util_set2(0, g.capacity() - 1, bitset::empt)
        , neighborhood(0, g.capacity() - 1, bitset::empt)
        , ordering_tmp(0, g.capacity() - 1, bitset::empt)
        , ordering_forbidden(0, g.capacity() - 1, bitset::empt)
        , ordering_removed_bs(0, g.capacity() - 1, bitset::empt)
        , expl_reps_extra(g.capacity())
        , altered_vertex_set(0, g.capacity() - 1, bitset::empt)
    {

        for (size_t i{0}; i < isrep.size(); ++i)
            isrep[i].initialise(0, g.capacity() - 1, gc::bitset::empt);

        stat.binds(this);
        ub = g.capacity()+1;

        expl_reps.clear();
        expl_reps.resize(g.capacity(), -1);

				varinfo.resize(s.nVars(), varinfo_t{-1,-1});

        for (int i = 0; i != g.capacity(); ++i) {
            if (!g.nodes.contain(i))
                continue;
            for (int j = 0; j != g.capacity(); ++j) {
                if (!g.nodes.contain(j))
                    continue;
                if (j < i) {
                    assert(vars[i][j] == vars[j][i]);
                    continue;
                }
                if (j == i)
                    continue;
                auto ijvar = vars[i][j];
                if (ijvar == minicsp::var_Undef)
                    continue;
                // if (varinfo.size() <= static_cast<size_t>(ijvar))
                //     varinfo.resize(ijvar + 1, varinfo_t{-1,-1});
                varinfo[ijvar] = {i, j};
                s.wake_on_lit(ijvar, this, nullptr);
                s.schedule_on_lit(ijvar, this);
            }
        }
        diffuv.initialise(0, g.capacity(), bitset::empt);
        diffvu.initialise(0, g.capacity(), bitset::empt);
        // cur_neighbors.initialise(0, g.capacity(), bitset::empt);

        if (opt.adaptive) {
            adaptive_callback
                = [this](const vec<Lit>&,
                      int) -> minicsp::Solver::clause_callback_result_t {
                this->run_expensive_bound = true;
                return minicsp::Solver::CCB_OK;
            };
            s.use_clause_callback(adaptive_callback);
        }

        if (opt.fillin == options::FILLIN_DECOMPOSE) {
            post_transitivity_decomposition();
        }

        // TODO
        // std::cout << "UB = " << ub << std::endl;
        // int k = ub-1;
        // int a = g.capacity() / k;
        //
        //
        //
        //                        mineq_constraint = new cons_pb(s, vars,
        // weights,
        // lb);
        //                        s.addConstraint(c);
        //
        // post_pb(Solver& s, std::vector<Var> const& vars,
        //               std::vector<int> const& weights, int lb);
        // int mineq = a * (g.capacity() - k * (a + 1) / 2);
        //                        mineq_constraint = new cons_pbvar(s, vbool, c,
        // rhs);
        //                        s.addConstraint(mineq_constraint);

        if (opt.prune)
            cf.switch_pruning();
        // cf.update_nn = opt.prune;

        if (opt.triangle_up)
            s.use_UP_callback([this](Lit p, const Clause& c) {
                return unit_propagate(p, c);
            });

        DO_OR_THROW(propagate(s));

    }

    std::ostream& print(Solver&, std::ostream& os) const
    {
        os << "coloring";
        return os;
    }

    std::ostream& printstate(Solver&, std::ostream& os) const
    {
        os << "coloring";
        return os;
    }

    void post_transitivity_decomposition()
    {
        int numclauses{0};
        // posts a clause assuming lit_Undef is false and ~lit_Undef is true
        auto post_safe = [&](std::vector<Lit>&& ps) {
            erase_if(ps, [](Lit l) { return l == lit_Undef; });
            for (Lit l : ps)
                if (l == ~lit_Undef)
                    return;
            s.addClause(std::move(ps));
            ++numclauses;
        };
        auto post_ternary_clauses = [&](int u, int v, int w) {
            Lit vu = Lit(vars[v][u]);
            Lit uw = Lit(vars[u][w]);
            Lit vw = Lit(vars[v][w]);
            post_safe({~Lit(vu), ~Lit(uw), Lit(vw)});
            post_safe({~Lit(vw), ~Lit(uw), Lit(vu)});
            post_safe({~Lit(vu), ~Lit(vw), Lit(uw)});
        };
        for (auto v : g.nodes) {
            util_set.copy(fg.origmatrix[v]);
            util_set.setminus_with(g.origmatrix[v]);
            for (auto u : util_set) {
                if (u < v)
                    continue;
                for (auto w : g.nodes) {
                    if (w < v)
                        continue;
                    if (fg.origmatrix[u].fast_contain(w)
                        && fg.origmatrix[v].fast_contain(w)) {
                        post_ternary_clauses(u, v, w);
                    }
                }
            }
        }
        std::cout << "[modeling] Decomposed transitivity constraint with "
                  << numclauses << " clauses\n";
    }

    minicsp::Solver::UP_result unit_propagate(Lit, const Clause& c)
    {
        sync_graph();
        int bestidx{-1};
        Lit other = c[0];
        auto s0 = sign(c[0]);
        varinfo_t info = varinfo[var(other)];
        int u = g.rep_of[info.u], v = g.rep_of[info.v];
        varinfo_t kinfo;
        int ku{-1}, kv{-1};
        for (int k = 2; k < c.size(); k++) {
            auto kvalue = s.value(c[k]);
            if (kvalue == l_False)
                continue;
            if (kvalue == l_True) {
                bestidx = k;
                break;
            }
            if (bestidx < 0) {
                bestidx = k;
            }
            if (sign(c[k]) != s0) {
                bestidx = k;
                break;
            }
            kinfo = varinfo[var(c[k])];
            ku = g.rep_of[kinfo.u];
            kv = g.rep_of[kinfo.v];
            if (ku != u || kv != v) {
                bestidx = k;
                break;
            }
        }
        if (ku == u && kv == v && s.value(c[bestidx]) != l_True
            && sign(c[bestidx]) == s0) {
#ifdef DEBUG_SUPERUP
            std::cout << "propportunity " << minicsp::print(s, &c) << "\n";
            std::cout << "unset literals: ";
            for (Lit l : c)
                if (s.value(l) != l_False)
                    std::cout << minicsp::lit_printer(s, l) << " ";
            std::cout << "\n";
            std::cout << "partition[" << v
                      << "] = " << print_container{g.partition[v]} << std::endl;
            std::cout << "partition[" << u
                      << "] = " << print_container{g.partition[u]} << "\n";
#endif
            util_set.clear();
            bool positive_clause{!sign(c[0])};
            reason.clear();
            for (Lit l : c) {
                if (s.value(l) == l_False)
                    reason.push(l);
                else {
                    kinfo = varinfo[var(l)];
                    util_set.fast_add(kinfo.u);
                    util_set.fast_add(kinfo.v);
                }
            }
            int urep{-1}, vrep{-1};
#ifdef DEBUG_SUPERUP
            std::cout << "vertices involved " << util_set << "\n";
            std::cout << "Reason: ";
#endif
            for (auto vp : g.partition[v]) {
                if (util_set.fast_contain(vp)) {
                    if (vrep < 0)
                        vrep = vp;
                    else {
#ifdef DEBUG_SUPERUP
                        std::cout << lit_printer(s, ~Lit(vars[vrep][vp]))
                                  << " ";
#endif
                        reason.push(~Lit(vars[vrep][vp]));
                    }
                }
            }
            for (auto up : g.partition[u]) {
                if (util_set.fast_contain(up)) {
                    if (urep < 0)
                        urep = up;
                    else {
#ifdef DEBUG_SUPERUP
                        std::cout << lit_printer(s, ~Lit(vars[urep][up]))
                                  << " ";
                        std::cout << std::flush;
#endif
                        reason.push(~Lit(vars[urep][up]));
                        assert(s.value(~Lit(vars[urep][up])) == l_False);
                    }
                }
            }
#ifdef DEBUG_SUPERUP
            std::cout << "\nEffect: " << lit_printer(s, Lit(vars[vrep][urep]))
                      << std::endl;
#endif
            if (s.value(vars[vrep][urep]) != l_Undef)
                return {bestidx};
            // assert(s.value(vars[vrep][urep]) == l_Undef);
            for (Lit l : reason)
                assert(s.value(l) == l_False);
            if (positive_clause)
                s.enqueueFill(Lit(vars[vrep][urep]), reason);
            else
                s.enqueueFill(~Lit(vars[vrep][urep]), reason);
        }
        return {bestidx};
    }

    // helper: because u merged with v and v merged with x, x
    // merged with u
    auto merge_3way(int u, int v, int x) -> Clause* {
        Var uxvar = vars[u][x];
        // if (opt.fillin && uxvar == var_Undef)
        //     return NO_REASON;
        if (s.value(uxvar) == l_True)
            return NO_REASON;
        if (u == v)
            return NO_REASON;
        if (!vars.contain(u, x))
            return NO_REASON;

        reason.clear();
        reason.push(~Lit(vars[u][v]));
        reason.push(~Lit(vars[v][x]));
        if (g.origmatrix[x].fast_contain(u)) {
            return s.addInactiveClause(reason);
        }
        DO_OR_RETURN(s.enqueueFill(Lit(uxvar), reason));
        return NO_REASON;
    };

    // helper: u is merged with v and v has an edge with x, so
    // u has an edge with x
    auto separate(int u, int v, int x) -> Clause* {
        Var uxvar = vars[u][x];
        if (opt.fillin && uxvar == var_Undef)
            return NO_REASON;
        if (g.origmatrix[u].fast_contain(x))
            return NO_REASON;
        if (u == v)
            return NO_REASON;
        if (!vars.contain(u, x))
            return NO_REASON;

        assert(u != x && v != x);
        reason.clear();
        reason.push(~Lit(vars[u][v]));
        if (!g.origmatrix[v].fast_contain(x)) {
            reason.push(Lit(vars[v][x]));
        }
        DO_OR_RETURN(s.enqueueFill(~Lit(vars[u][x]), reason));
        return NO_REASON;
    };

    Clause* transitive_separate(
        int u, int v, const bitset& diffvu, const bitset& vpart)
    {
        std::cout << "transitive separate " << u << "--" << v << "\n";
        std::cout << "partition[" << u
                  << "] = " << print_container<std::vector<int>>{g.partition[u]} << "\n";
        std::cout << "partition[" << v
                  << "] = " << print_container<std::vector<int>>{g.partition[v]} << "\n";
        std::cout << u << " will gain neighbors " << diffvu << "\n";
        needbfs.clear();
        // for each partition to which we must add edges, first find
        // at least one edge that can be fixed immediately, which will
        // serve as the root of the BFS. This will actually set all
        // such edges and keep one of them as the root
        for (auto up : g.partition[u]) {
            util_set.copy(diffvu);
            util_set.intersect_with(fg.origmatrix[up]);
            for (auto n : util_set) {
                if (fg.origmatrix[up].fast_contain(v)
                    && fg.origmatrix[n].fast_contain(v)) {
                    std::cout << "Immediate separation " << up << "--" << n
                              << " via rep = " << v << "\n";
                    DO_OR_RETURN(separate(up, v, n));
                    bfsroots[n] = {up, n};
                } else {
                    util_set2.copy(vpart);
                    util_set2.intersect_with(fg.origmatrix[up]);
                    util_set2.intersect_with(fg.origmatrix[n]);
                    if (util_set2.empty()) {
                        needbfs.fast_add(g.rep_of[n]);
                        // bfsroots[n] = {up, n};
                        std::cout << "Will set " << up << "--" << n
                                  << " using BFS\n";
                        continue;
                    }
                    int via = util_set2.min();
                    std::cout << "Immediate step of separation " << up << "--"
                              << n << " via " << via << "\n";
                    DO_OR_RETURN(separate(up, via, n));
                    bfsroots[n] = {up, n};
                    std::cout << "\t\tWill start BFS from " << up
                              << " to reach " << n << "\n";
                }
            }
        }

        for (auto v : needbfs) {
            auto e = bfsroots[v];
            std::cout << "\n****\tBFSing " << e.first << "-" << e.second
                      << "\n";
            bfs.bfs(e.second, e.first, [&](int x, int w) -> Clause* {
                std::cout << "transitive_separate: checking edge " << w << "--"
                          << x << "\n";
                bool xinv = bfs.vpart.fast_contain(x);
                bool winv = bfs.vpart.fast_contain(w);
                std::cout << "inv: " << winv << xinv << "\n";
                // std::cout << "vpart = " << bfs.vpart << "\n";
                // std::cout << "upart = " << bfs.upart << "\n";
                if (xinv == winv)
                    return NO_REASON;
                if (xinv) {
                    std::cout << "\tseparate(" << x << "," << bfs.vparent[w]
                              << "," << w << ")\n";
                    DO_OR_RETURN(separate(x, bfs.vparent[w], w));
                } else {
                    std::cout << "here\n";
                    std::cout << "\tseparate(" << x << "," << bfs.uparent[w]
                              << "," << w << ")\n";
                    DO_OR_RETURN(separate(x, bfs.uparent[w], w));
                }
                return NO_REASON;
            });
        }

        return NO_REASON;
    }

    Clause *merge_fillin(Lit l)
    {
        auto info = varinfo[var(l)];
        auto u{g.rep_of[info.u]}, v{g.rep_of[info.v]};
        if (u == v)
            return NO_REASON;

        if (opt.fillin == options::FILLIN_DECOMPOSE) {
            g.merge(u, v);
#ifdef UPDATE_FG
            fg.merge(u, v);
#endif
            return NO_REASON;
        }

        std::cout << "\n\nmerge " << info.u << "--" << info.v << "\n";
        std::cout << "partition[" << u
                  << "] = " << print_container<std::vector<int>>{g.partition[u]} << "\n";
        std::cout << "partition[" << v
                  << "] = " << print_container<std::vector<int>>{g.partition[v]} << "\n";

        // same as in plain version of wake()
        diffuv.copy(g.matrix[u]);
        diffuv.setminus_with(g.matrix[v]);
        diffvu.copy(g.matrix[v]);
        diffvu.setminus_with(g.matrix[u]);

        // set to true all the variables between the two partitions
        bfs.bfs(info.u, info.v, [&](int w, int x) -> Clause* {
            if ((x == info.u || x == info.v) && (w == info.u || w == info.v))
                return NO_REASON;
            // followed edge from w to x
            std::cout << "merge: checking edge " << w << "->" << x << "\n";
            bool xinv = bfs.vpart.fast_contain(x);
            bool winv = bfs.vpart.fast_contain(w);
            std::cout << "inv: " << winv << xinv << "\n";
            std::cout << "vpart = " << bfs.vpart << "\n";
            std::cout << "upart = " << bfs.upart << "\n";
            if (xinv == winv)
                return NO_REASON;
            if (xinv) {
                std::cout << "\t1: merge(" << x << "," << bfs.vparent[w] << ","
                          << w << ")\n";
                assert(bfs.vparent[w] != x || x == bfs.rootv);
                if (x == bfs.rootv)
                    DO_OR_RETURN(merge_3way(x, bfs.rootu, w));
                else
                    DO_OR_RETURN(merge_3way(x, bfs.vparent[w], w));
            } else {
                std::cout << "\t2: merge(" << x << "," << bfs.uparent[w] << ","
                          << w << ")\n";
                assert(bfs.uparent[w] != x || x == bfs.rootu);
                if (x == bfs.rootu)
                    DO_OR_RETURN(merge_3way(x, bfs.rootv, w));
                else
                    DO_OR_RETURN(merge_3way(x, bfs.uparent[w], w));
            }
            return NO_REASON;
        });

        // now use the neighborhood differences to separate each
        // vertex from new neighbors.
        vpartcopy.copy(bfs.vpart); // this will be destroyed by the
                                   // bfs calls in transitive_separate()
        DO_OR_RETURN(transitive_separate(v, info.u, diffuv, bfs.upart));
        DO_OR_RETURN(transitive_separate(u, info.v, diffvu, vpartcopy));

        g.merge(u, v);
#ifdef UPDATE_FG
        fg.merge(u, v);
#endif
        return NO_REASON;
    }

    Clause *separate_fillin(Lit l)
    {
        auto info = varinfo[var(l)];
        auto u{g.rep_of[info.u]}, v{g.rep_of[info.v]};
        if (g.matrix[u].fast_contain(v))
            return NO_REASON;

        if (opt.fillin == options::FILLIN_DECOMPOSE) {
            g.separate(u, v);
#ifdef UPDATE_FG
            fg.separate(u, v);
#endif
            return NO_REASON;
        }

        bfs.bfs(info.u, info.v, [&](int x, int w) -> Clause* {
            std::cout << "separate: checking edge " << w << "--" << x << "\n";
            if ((x == info.u || x == info.v) && (w == info.u || w == info.v))
                return NO_REASON;
            std::cout << "pass 1\n";
            bool xinv = bfs.vpart.fast_contain(x);
            bool winv = bfs.vpart.fast_contain(w);
            std::cout << "inv: " << winv << xinv << "\n";
            std::cout << "vpart = " << bfs.vpart << "\n";
            std::cout << "upart = " << bfs.upart << "\n";
            if (xinv == winv)
                return NO_REASON;
            if (xinv) {
                std::cout << "vparent[" << w << "] = " << bfs.vparent[w]
                          << "\n";
                DO_OR_RETURN(separate(x, bfs.uparent[x], w));
            } else {
                std::cout << "uparent[" << w << "] = " << bfs.uparent[w]
                          << "\n";
                DO_OR_RETURN(separate(x, bfs.vparent[x], w));
            }
            return NO_REASON;
        });
        g.separate(u, v);
#ifdef UPDATE_FG
        fg.separate(u, v);
#endif
        return NO_REASON;
    }

    Clause *wake_fillin(Lit l)
    {
        if (!sign(l))
            return merge_fillin(l);
        else
            return separate_fillin(l);
    }

    Clause* wake(Solver& s, Lit l)
    {
        sync_graph();

        if (opt.cliquealg == options::INCREMENTAL) {
            auto info = varinfo[var(l)];

            // std::cout << " " << info.u << "|" << info.v;

            if (opt.prune and g.nodeset.fast_contain(info.u)
                and g.nodeset.fast_contain(info.v)) {
                arc e{info.u, info.v};
                changed_edges.push_back(e);
                // changed_edges.push_back(info.v);
            }

            // std::cout << "add";
            if (!altered_vertex_set.fast_contain(info.u)) {
                altered_vertices.push_back(info.u);
                altered_vertex_set.fast_add(info.u);

                // if(g.nodeset.fast_contain(info.u))
                //      std::cout << " " << info.u;
            }
            if (!altered_vertex_set.fast_contain(info.v)) {
                altered_vertices.push_back(info.v);
                altered_vertex_set.fast_add(info.v);

                // if(g.nodeset.fast_contain(info.v))
                //      std::cout << " " << info.v;
            }
            // std::cout << std::endl;
        }

        if (opt.fillin)
            return wake_fillin(l);

        auto info = varinfo[var(l)];
        auto u{g.rep_of[info.u]}, v{g.rep_of[info.v]};

        if (!sign(l)) {
            // merging u and v

            if (u == v) {
                // already merged
                return NO_REASON;
            }

            // before merging, keep track of the differences in neighborhood
            diffuv.copy(g.matrix[u]);
            diffuv.setminus_with(g.matrix[v]);
            diffvu.copy(g.matrix[v]);
            diffvu.setminus_with(g.matrix[u]);

            // merge each partition with the other partition's
            // representative
            for (auto vp : g.partition[v]) {
                if (vp == info.v)
                    continue;
                DO_OR_RETURN(merge_3way(info.u, info.v, vp));
            }

            for (auto up : g.partition[u]) {
                if (up == info.u)
                    continue;
                DO_OR_RETURN(merge_3way(info.v, info.u, up));
            }

            // then merge all the rest, using the fact that each is
            // merged with the representative
            for (auto vp : g.partition[v]) {
                if (vp == info.v)
                    continue;

                for (auto up : g.partition[u]) {
                    if (up == info.u)
                        continue;
                    DO_OR_RETURN(merge_3way(vp, info.u, up));
                }
            }

            // now separate each vertex from potential new
            // neighbors. Here, we only use the original neighborhood
            // of the representative of each partition to determine
            // whether to generate a binary or ternary clause. We
            // could try to find a vertex in the partition that was
            // originally adjacent to the new neighbor, but it's too
            // much work (linear, even if it's mediated by a bitset)

            for (auto vp : g.partition[v]) {
                for (auto n : diffuv)
                    DO_OR_RETURN(separate(vp, info.u, n));
            }
            for (auto up : g.partition[u]) {
                for (auto n : diffvu)
                    DO_OR_RETURN(separate(up, info.v, n));
            }

            g.merge(u, v);
            dynamic_reps[u] = info.u;
            dynamic_reps[v] = info.v; // the partition is not there
                                      // anymore, but we may use this
                                      // when we backtrack
        } else {

            if (g.matrix[u].fast_contain(v))
                return NO_REASON;

            // separate each vertex in each partition from the
            // representative of the other partition
            for (auto up : g.partition[u])
                DO_OR_RETURN(separate(up, info.u, info.v));
            for (auto vp : g.partition[v])
                DO_OR_RETURN(separate(vp, info.v, info.u));

            // and now all-to-all
            for (auto up : g.partition[u])
                for (auto vp : g.partition[v])
                    DO_OR_RETURN(separate(up, info.u, vp));

            g.separate(u, v);
        }

        return NO_REASON;
    }

    // If opt.learning is CHOOSE_POSITIVE, itfinds the variable
    // between two partitions with max VSIDS and sets expl_reps of the
    // two partitions so that they choose this variable. Otherwise,
    // just sets expl_reps to be equal to the representative of the
    // partition (static or dynamic)
    //
    // If expl_reps for one of the two partitions is already set, its
    // choice is consistent with that value.
    //
    // returns true if it found an edge that was there in the original
    // graph (which is the best match there can be, because we do not
    // need to put a literal for this edge in the clause)
    bool bestmatch(int u, int v, const std::vector<int>& ubag,
        const std::vector<int>& vbag, int up, int& maxvar, double& maxactivity)
    {
        assert(u == g.rep_of[u]);
        assert(v == g.rep_of[v]);
        if (opt.learning != options::CHOOSE_POSITIVE) {
            if (opt.learning == options::DYNAMIC_REPS) {
                int ur = dynamic_reps[u];
                int vr = dynamic_reps[v];
                // use the dynamic representative, except if it's not
                // in the same partition anymore (can happen after
                // backtracking)
                if (g.rep_of[ur] != u)
                    ur = g.rep_of[u];
                if (g.rep_of[vr] != v)
                    vr = g.rep_of[v];
                expl_reps[u] = ur;
                expl_reps[v] = vr;
            } else {
                expl_reps[u] = g.rep_of[u];
                expl_reps[v] = g.rep_of[v];
            }
            if (g.origmatrix[expl_reps[u]].fast_contain(expl_reps[v]))
                maxvar = var_Undef;
            else
                maxvar = vars[expl_reps[u]][expl_reps[v]];
            return true;
        }
        int bestu{-1}, bestv{-1};
        for (auto vp : vbag) {
            if (g.origmatrix[up].fast_contain(vp)) {
                expl_reps[u] = up;
                expl_reps[v] = vp;
                maxvar = var_Undef;
                return true;
            } else {
                auto var = vars[up][vp];
                auto varact = s.var_activity(var);
                assert(var != var_Undef);
                assert(s.value(var) == l_False);
                if (varact > maxactivity) {
                    maxactivity = varact;
                    maxvar = var;
                    bestu = up;
                    bestv = vp;
                }
            }
        }
        if (bestu >= 0) {
            expl_reps[u] = bestu;
            expl_reps[v] = bestv;
        }
        return false;
    }

    // helper for explain_positive_clique and for myc
    // explanations. Explains why an edge is there between two bags,
    // taking into account chosen (or not-yet-chosen) reps and fillin
    void explain_edge(int u, int v)
    {
        assert(g.rep_of[u] == u);
        assert(g.rep_of[v] == v);
        assert(g.matrix[u].fast_contain(v));
        assert(g.matrix[v].fast_contain(u));
        auto ubag = &g.partition[u];
        auto vbag = &g.partition[v];
        double maxactivity = -1.0;
        Var maxvar = var_Undef;
        if (expl_reps[v] >= 0 && expl_reps[u] < 0) {
            using std::swap;
            swap(u, v);
            swap(ubag, vbag);
        }
        if (expl_reps[v] < 0 && expl_reps[u] >= 0) {
            // find a rep for vbag only
            auto up = expl_reps[u];
            bestmatch(u, v, *ubag, *vbag, up, maxvar, maxactivity);
        } else if (expl_reps[v] < 0 && expl_reps[u] < 0) {
            // find a rep for both vbag and ubag
            for (auto up : *ubag) {
                if (bestmatch(u, v, *ubag, *vbag, up, maxvar, maxactivity))
                    break;
            }
        } else {
            auto ur = expl_reps[u];
            auto vr = expl_reps[v];
            maxvar = vars[ur][vr];
            assert(maxvar == var_Undef || s.value(maxvar) == l_False);
        }
        bool found = maxvar != var_Undef
            || g.origmatrix[expl_reps[u]].fast_contain(expl_reps[v]);
        assert(found || opt.fillin);
        if (maxvar != var_Undef) {
            assert(s.value(maxvar) == l_False);
            reason.push(Lit(maxvar));
            return;
        }
        if (found)
            return;

        assert(opt.fillin);

        auto find_variable
            = [&](auto& vertices, auto& N) -> std::tuple<bool, int, int> {
            for (auto up : vertices) {
                util_set2.copy(N);
                util_set2.intersect_with(fg.origmatrix[up]);
                if (!util_set2.empty()) {
                    int vp = util_set2.min();
                    maxvar = vars[up][vp];
                    CCOUT << "Using edge " << vp << "--" << up << "\n"
                          << "expl_reps[" << u << "] = " << expl_reps[u] << "\n"
                          << "expl_reps_extra[" << u
                          << "] = " << print_container<std::vector<int>>{expl_reps_extra[u]}
                          << "\n"
                          << "expl_reps[" << v << "] = " << expl_reps[v] << "\n"
                          << "expl_reps_extra[" << v
                          << "] = " << print_container<std::vector<int>>{expl_reps_extra[v]}
                          << "\n";
                    // if edge is in original graph, no var
                    if (maxvar != var_Undef) {
                        assert(s.value(Lit(maxvar)) == l_False);
                        reason.push(Lit(maxvar));
                    }
                    return {true, up, vp};
                }
            }
            return {false, var_Undef, var_Undef};
        };

        auto use_extra_rep = [&](int v, int v2) {
            CCOUT << "Using extra rep " << v2 << " for " << v << "\n";
            if (expl_reps[v] == -1) {
                expl_reps[v] = v2;
                return;
            }
            int v1 = expl_reps[v];
            // XXX: this is a linear search, which could
            // be expensive if we somehow have to use many
            // extra reps
            if (v1 == v2
                || find(expl_reps_extra[v1].begin(), expl_reps_extra[v1].end(),
                       v2)
                    != expl_reps_extra[v1].end())
                return;
            // fuck, we have to explain why v2 is in the same
            // partition as expl_reps[v]
            if (fg.origmatrix[v1].fast_contain(v2)) {
                reason.push(~Lit(vars[v1][v2]));
                assert(vars[v1][v2] != var_Undef);
                assert(s.value(~Lit(vars[v1][v2])) == l_False);
                expl_reps_extra[v].push_back(v2);
            } else {
                // breadth first search
                found = false;
                CCOUT << "Trying to explain that " << v1 << " and " << v2
                      << " are in the same partition\n";
                bfs.bfs(v1, [&](int w1, int w2) {
                    if (w2 == v2) {
                        CCOUT << "Arrived!\n";
                        found = true;
                        do {
                            CCOUT << "using edge " << w2 << "--"
                                  << bfs.parent[w2] << "\n";
                            Var var = vars[w2][bfs.parent[w2]];
                            assert(var != var_Undef);
                            assert(s.value(var) == l_True);
                            reason.push(~Lit(var));
                            CCOUT << "\tusing " << lit_printer(s, ~Lit(var))
                                  << "\n";
                            w2 = bfs.parent[w2];
                        } while (w2 != v1);
                        return true;
                    }
                    return false;
                });
                if (!found)
                    std::cout << std::endl;
                assert(found);
            }
            CCOUT << "expl_reps[" << v << "] = " << expl_reps[v] << "\n"
                  << "expl_reps_extra[" << u
                  << "] = " << print_container<std::vector<int>>{expl_reps_extra[v]} << "\n";
        };

        CCOUT << "explaining edge " << u << "--" << v << " the hard way\n"
              << "partitions[u] = " << print_container<std::vector<int>>{g.partition[u]} << "\n"
              << "partitions[v] = " << print_container<std::vector<int>>{g.partition[v]} << "\n";

        // first, try to find a variable
        // between the existing extra reps
        util_set.clear();
        for (auto vp : expl_reps_extra[v])
            util_set.fast_add(vp);
        std::tuple<bool, int, int> rv{false, -1, -1};
        if (!util_set.empty())
            rv = find_variable(expl_reps_extra[u], util_set);

        using std::get;
        if (get<0>(rv))
            return;
        // now try to find a variable between one of the
        // ureps and the rest of the partition of v
        for (auto vp : *vbag)
            util_set.fast_add(vp);
        rv = find_variable(expl_reps_extra[u], util_set);
        if (get<0>(rv)) {
            use_extra_rep(v, get<2>(rv));
        } else {
            // worst case: try to find a variable between any
            // of the vertices in vbag and any of the vertices
            // in ubag
            rv = find_variable(*ubag, util_set);
            assert(get<0>(rv));
            use_extra_rep(u, get<1>(rv));
            use_extra_rep(v, get<2>(rv));
        }
        // std::cout << "current reason " << minicsp::print(s, &reason)
        //           << "\n";
    }

    // fills in reason with explanation for clq. Also generates
    // this->culprit and fills in this->expl_reps, unless chosen_reps
    // is true, in which case it assumes those are already set
    // externally. In any case, it might also populate
    // expl_reps_extra, if opt.fillin in not NONE
    void explain_positive_clique(const bitset& clq, bool chosen_reps)
    {
        culprit.clear();
        std::copy(begin(clq), end(clq), back_inserter(culprit));
        if (opt.learning == options::CHOOSE_POSITIVE) {
            // sort by increasing partition size
            std::sort(begin(culprit), end(culprit), [&](auto u, auto v) {
                return g.partition[u].size() < g.partition[v].size();
            });
        }

        assert(!culprit.empty());
        CCOUT << "Explaining " << print_container<std::vector<int>>{culprit} << "\n";

        if (!chosen_reps) {
            expl_reps.clear();
            expl_reps.resize(g.capacity(), -1);
        }
        for (auto v : culprit)
            expl_reps_extra[v].clear();

        for (size_t i = 0; i != culprit.size() - 1; ++i)
            for (size_t j = i + 1; j != culprit.size(); ++j) {
                auto u = culprit[i], v = culprit[j];
                explain_edge(u,v);
            }
        CCOUT << "Finished with reason " << minicsp::print(s, &reason) << "\n";
        for (Lit l : reason) {
            assert(var(l) != var_Undef);
            assert(s.value(l) == l_False);
        }
    }

    Clause* explain_positive()
    {
        auto maxidx{std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques))};

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            maxidx = mf.explanation_clique;
        }

        // explain the base clique
        reason.clear();
        explain_positive_clique(cf.cliques[maxidx], false);

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            // for (auto v : mf.explanation_subgraph.nodes) {
            for (auto i{culprit.size()};
                 i < mf.explanation_subgraph.nodes.size(); ++i) {
                auto v{mf.explanation_subgraph.nodes[i]};
                neighborhood.copy(mf.explanation_subgraph.matrix[v]);
                for (auto u : neighborhood) {
                    if (mf.explanation_subgraph.nodes.index(v)
                        > mf.explanation_subgraph.nodes.index(u)) {
                        explain_edge(v, u);
                    }
                }
            }
        }

        return s.addInactiveClause(reason);
    }


    Clause* explain()
    {
        switch (opt.learning) {
        case options::NO_LEARNING:
            return INVALID_CLAUSE;
        case options::NAIVE_POSITIVE:
        case options::CHOOSE_POSITIVE:
        case options::DYNAMIC_REPS:
            return explain_positive();
        default:
            assert(0);
            return INVALID_CLAUSE;
        }
    }

    Clause* contractWhenNIncluded()
    {
        bool some_propagation = true;
        while (some_propagation) {
            some_propagation = false;
            for (auto u : g.nodes) {
                for (auto v : g.nodes) {
                    if (u != v && !g.origmatrix[u].fast_contain(v)
                        && s.value(vars[u][v]) == l_Undef) {
                        neighborhood.copy(g.matrix[v]);
                        neighborhood.setminus_with(g.matrix[u]);
                        if (!neighborhood.intersect(g.nodeset)) {
                            // N(v) <= N(U)s
                            some_propagation = true;
                            ++stat.num_neighborhood_contractions;

                            // explanation: the edges (ON(v) \ ON(u)) x u
                            neighborhood.copy(g.origmatrix[v]);
                            neighborhood.setminus_with(g.origmatrix[u]);
                            neighborhood.intersect_with(g.nodeset);

                            reason.clear();
                            for (auto w : neighborhood) {
                                reason.push(Lit(vars[u][w]));
                            }
                            DO_OR_RETURN(
                                s.enqueueFill(Lit(vars[u][v]), reason));
                        }
                    }
                }
            }
        }
        return NO_REASON;
    }

    Clause* lastColorPruning()
    {

#ifdef _DEBUG_PRUNING
        cf.notify_ub(ub);
#endif

#ifdef _FULL_PRUNING
        erase_if(changed_edges, [&](arc e) {
            return !g.nodeset.fast_contain(e[0])
                or !g.nodeset.fast_contain(e[1]);
        });

        assert(changed_vertices.size() == 0);
        // changed_vertices.clear();
        altered_vertex_set.clear();

        for (auto e : changed_edges) {
            for (int i = 0; i < 2; ++i)
                if (!altered_vertex_set.fast_contain(e[i])) {
                    changed_vertices.push_back(e[i]);
                    altered_vertex_set.fast_add(e[i]);
                }
        }
#endif

        cf.compute_pruning(ub - 1, changed_edges, changed_vertices);

				// std::cout << opt.enurp << std::endl;


        if (opt.enurp and cf.new_tight_cliques.size() > 0) {
					
            for (int i = 0; i < s.nVars(); ++i) {
							
							// std::cout << i << std::endl;
                auto info{varinfo[i]};
                if (!info.empty() and s.value(i) == l_Undef) {
                    auto u{g.rep_of[info.u]};
                    auto v{g.rep_of[info.v]};
                    if (!g.origmatrix[v].fast_contain(u)) {

                        for (auto i : cf.new_tight_cliques) {
                            if (cf.cliques[i].fast_contain(u)
                                or cf.cliques[i].fast_contain(v)
                                or (cf.neighbor_count[i][u]
                                       + cf.neighbor_count[i][v])
                                    < cf.clique_sz[i])
                                continue;

                            reason.clear();

                            bool included{true};
                            for (auto w : cf.cliques[i]) {
                                if (g.matrix[v].fast_contain(w)) {
                                    explain_edge(w, v);
                                } else if (g.matrix[u].fast_contain(w)) {
                                    explain_edge(u, w);
                                } else {
                                    included = false;
                                    break;
                                }
                            }

                            if (included) {
                                explain_positive_clique(cf.cliques[i], false);
                                DO_OR_RETURN(
                                    s.enqueueFill(~Lit(vars[u][v]), reason));
                            }
                        }
                    }
                }
            }
        }
				cf.new_tight_cliques.clear();

#ifdef _FULL_PRUNING
        changed_edges.clear();
        changed_vertices.clear();
#endif


#ifdef _DEBUG_PRUNING
        for (auto v : g.nodes) {
            for (auto i{0}; i < cf.num_cliques; ++i) {
                assert(cf.num_neighbors_of[i][v] == cf.neighbor_count[i][v]);
            }
        }
        assert(cf.pruning.size() == cf.col_pruning.size());
#endif

        while (cf.col_pruning.size() > 0) {
            auto u = cf.col_pruning.back();
            cf.col_pruning.pop_back();

            auto i = cf.col_pruning.back();
            cf.col_pruning.pop_back();

            util_set.copy(cf.cliques[i]);
            util_set.setminus_with(g.matrix[u]);

            auto w{util_set.min()};
            assert(util_set.size() == 1);

            reason.clear();
            explain_positive_clique(cf.cliques[i], false);
            for (auto v1 : cf.cliques[i]) {
                if (v1 != w)
                    explain_edge(v1, u);
            }
            DO_OR_RETURN(s.enqueueFill(Lit(vars[u][w]), reason));
        }

        return NO_REASON;
    }

    // find an upper bound u for the maximum IS. n/ceil(u) is a lower
    // bound for the chromatic number
    Clause* propagate_is()
    {
        create_ordering();
        ccf.find_clique_cover(heuristic);
        if (std::ceil(g.nodes.size() / static_cast<double>(ccf.num_cliques))
            >= ub) {
            std::cout << "\tIS ub = " << ccf.num_cliques
                      << " |V| = " << g.nodes.size() << " coloring lb = "
                      << std::ceil(g.nodes.size()
                             / static_cast<double>(ccf.num_cliques))
                      << "\n";
            reason.clear();
            for (int i = 0; i != ccf.num_cliques; ++i)
                explain_positive_clique(ccf.cliques[i], false);
            return s.addInactiveClause(reason);
        }
        return NO_REASON;
    }

    void create_ordering() { create_ordering(g.nodes); }

    template <typename NodeVec> void create_ordering(const NodeVec& hnodes)
    {

        // recompute the degenracy order
        if (opt.ordering == options::DYNAMIC_DEGENERACY) {
            assert(false);
        } else if (opt.ordering == options::DEGENERACY
            or opt.ordering == options::INVERSE_DEGENERACY) {
            assert(false);
        } else if (opt.ordering == options::PARTITION) {

            if (opt.ordering_low_degree == options::PREPROCESSING_ORDERING) {

                bool removed_some{false};
                ordering_removed_bs.clear();
                ordering_removed.clear();
                do {

                    removed_some = false;
                    ordering_forbidden.clear();
                    for (auto v : hnodes) {

                        std::cout << "failed assert after checking " << v
                                  << std::endl;

                        if (ordering_forbidden.fast_contain(v)
                            || ordering_removed_bs.fast_contain(v))
                            continue;
                        ordering_tmp.copy(g.matrix[v]);
                        ordering_tmp.intersect_with(g.nodeset);
                        if (ordering_tmp.size()
                            >= static_cast<size_t>(std::max(bestlb, *lastlb)))
                            continue;
                        removed_some = true;
                        ordering_removed_bs.fast_add(v);
                        ordering_removed.push_back(v);
                        ordering_forbidden.union_with(g.matrix[v]);
                        ordering_forbidden.fast_add(v);
                    }
                } while (removed_some);
            } else if (opt.ordering_low_degree == options::DEGREE_ORDERING) {
                ordering_removed_bs.clear();
                ordering_removed.clear();
                for (auto v : hnodes) {
                    ordering_tmp.copy(g.matrix[v]);
                    ordering_tmp.intersect_with(g.nodeset);
                    if (ordering_tmp.size()
                        >= static_cast<size_t>(std::max(bestlb, *lastlb)))
                        continue;
                    ordering_removed_bs.fast_add(v);
                    ordering_removed.push_back(v);
                }
            }

            // sort by partition size
            heuristic.clear();
            for (auto v : hnodes)
                if (!ordering_removed_bs.fast_contain(v))
                    heuristic.push_back(v);

            std::sort(heuristic.begin(), heuristic.end(),
                [&](const int x, const int y) {
                    return (g.partition[x].size() > g.partition[y].size()
                        || (g.partition[x].size() == g.partition[y].size()
                               && x < y));
                });

            std::reverse(begin(ordering_removed), end(ordering_removed));
            for (auto v : ordering_removed)
                heuristic.push_back(v);

            assert(heuristic.size() == hnodes.size());
        } else {

            // no ordering
            heuristic.clear();
            std::copy(
                hnodes.begin(), hnodes.end(), std::back_inserter(heuristic));
        }
    }

    Clause* propagate_IS_constraints()
    {
        for (size_t i = 0; i < isconses.size(); ++i) {
            isrep[i].clear();
            for (auto v : isconses[i].vs) {
                isrep[i].fast_add(g.rep_of[v]);
            }
        }

        for (int i = 0; i != cf.num_cliques; ++i) {
            if (cf.clique_sz[i] < ub - 1)
                continue;

            for (size_t j = 0; j < isconses.size(); ++j) {
                const auto& c{isconses[j]};
                const auto& r{isrep[j]};
                util_set.copy(r);
                util_set.intersect_with(cf.cliques[i]);
                if (static_cast<int>(util_set.size()) >= ub - 1) {
                    // we need to choose the vertices of the
                    // constraint as representatives of the
                    // partitions participating in the clique
                    expl_reps.clear();
                    expl_reps.resize(g.capacity(), -1);
                    for (auto v : c.vs) {
                        if (!util_set.fast_contain(g.rep_of[v]))
                            continue;
                        expl_reps[g.rep_of[v]] = v;
                    }

                    reason.clear();
                    explain_positive_clique(util_set, true);
                    return s.addInactiveClause(reason);
                }
            }
        }
        return NO_REASON;
    }

    Clause* propagate(Solver&) final
    {
        // check_consistency();

        // std::cout << "\nPROPAGATE (" << s.decisionLevel() << ")\n";

        int lb{0};

        bound_source = options::CLIQUES;

        create_ordering();

        if (opt.cliquealg == options::ALL_CLIQUES || !saved_cliques) {
            lb = cf.find_cliques(heuristic, lb, ub, opt.cliquelimit);
            if (opt.cliquealg == options::INCREMENTAL) {
                cf.filter_cliques(lb, ub - 1);
                saved_cliques = true;
            }
        } else {
            erase_if(altered_vertices,
                [&](int v) { return !g.nodeset.fast_contain(v); });
            create_ordering(altered_vertices);
            lb = cf.extend_cliques(heuristic, lb, ub);
            cf.filter_cliques(lb - opt.cliquemargin, ub - 1);
            altered_vertices.clear();
            altered_vertex_set.clear();
        }

        stat.notify_nclique(cf.num_cliques);

        auto clique_bound{true};
        if (lb < ub
            && (s.decisionLevel() == 0 || !opt.adaptive
                   || run_expensive_bound)) {
            run_expensive_bound = false;
            bound_source = opt.boundalg;
            auto mlb{lb};

            // auto better{lb < mlb};

            if (opt.boundalg == options::FULLMYCIELSKI) {
                mlb = mf.full_myciel(lb, ub, s, vars);
            } else if (opt.boundalg == options::MAXMYCIELSKI) {
                mlb = mf.improve_cliques_larger_than(lb, lb, ub, s, vars);
            } else if (opt.boundalg == options::GREEDYMYCIELSKI) {
                mlb = mf.improve_greedy(lb - 1, lb, ub, s, vars);
            }

            clique_bound = (lb == mlb);

            stat.notify_bound_delta(lb, mlb);

            // if(stat.num_bound_delta%100 == 0)
            //              stat.display(std::cout);

            lb = mlb;
        }

        // if(lb == ub-1)
        //      std::cout << "prune? (" << s.decisionLevel() << ")\n";

        // std::cout << "LOWERBOUND (" << lb << " / " << bestlb << ")\n";

        *lastlb = lb;

        bool use_global_bound = false;
        if (lb <= bestlb) {
            lb = bestlb;
            use_global_bound = true;
        }

        if (s.decisionLevel() == 0 and lb > bestlb) {
            if (lb <= ub) {
                // this may not be the case when the instance is UNSAT
                bestlb = lb;
                stat.display(std::cout);
            }
        }

        if (cf.num_cliques == 1 && opt.cliquealg == options::ALL_CLIQUES)
            assert(g.nodes.size() == cf.cliques[0].size());

        if (lb >= ub) {
            if (use_global_bound) {
                reason.clear();
                return s.addInactiveClause(reason);
            }

            // std::cout << "FAIL! (" << s.decisionLevel() << ")\n";

            return explain();
        } else if (clique_bound and opt.prune) {
            if (lb == ub - 1)
                DO_OR_RETURN(lastColorPruning());
            else
                changed_edges.clear();
        }

        if (opt.indset_lb)
            DO_OR_RETURN(propagate_is());

        // check local constraints
#ifdef DEBUG_IS
        std::cout << "\npropag: " << lb << " " << ub
                  << " dlvl = " << s.decisionLevel() << std::endl;
#endif

        if (opt.indset_constraints and lb >= ub - 1 and lb >= ub - 1) {
            DO_OR_RETURN(propagate_IS_constraints());
        }

        // auto sol{gc::brelaz_color(g,true)};
        // int ncol{*max_element(begin(sol), end(sol)) + 1};
        // if(ncol == lb){
        //      std::cout << "STOP! (" << s.decisionLevel() << ")\n";
        // }

        return NO_REASON;
    }

    void check_consistency()
    {
        g.check_consistency();

        for (auto v : g.nodes) {
            assert(!g.partition[v].empty());
            assert(g.partition[v][0] == v);
            bitset bs(g.matrix[v]);
            for (auto u : g.partition[v]) {
                assert(g.rep_of[u] == v);
                bs.setminus_with(g.origmatrix[u]);
            }
            for (auto u : bs) {
                if (opt.fillin && !fg.origmatrix[v].fast_contain(u))
                    continue;
                assert(vars[v][u] != var_Undef);
                if (s.value(vars[v][u]) != l_False)
                    std::cout << "Partition " << v
                              << " = "
                              // << print_container(g.partition[v])
                              << " has extra neighbor " << u << std::endl;
            }
            for (auto u : bs) {
                if (opt.fillin && !fg.origmatrix[v].fast_contain(u))
                    continue;
                assert(s.value(vars[v][u]) == l_False);
            }
        }

        bool failed{false};
        for (auto& vec : vars.vars)
            for (auto p : vec.second) {
                auto x = p.second;
                if (x == var_Undef)
                    continue;
                auto info = varinfo[x];
                if (info.u < 0 || info.u == info.v)
                    continue;
                if (s.value(x) == l_False) {
                    if (!g.matrix[g.rep_of[info.u]].contain(g.rep_of[info.v])) {
                        std::cout << info.u << "-" << info.v
                                  << " not connected, var is false\n";
                        failed = true;
                    }
                    if (!g.matrix[g.rep_of[info.v]].contain(g.rep_of[info.u])) {
                        std::cout << info.v << "-" << info.u
                                  << " not connected, var is false\n";
                        failed = true;
                    }
                }
                if (s.value(x) == l_True) {
                    if (g.rep_of[info.u] != g.rep_of[info.u]) {
                        std::cout << info.u << "-" << info.v
                                  << " not merged, var is true\n";
                        failed = true;
                    }
                }
                if (s.value(x) == l_Undef) {
                    int ur = g.rep_of[info.u];
                    int vr = g.rep_of[info.v];
                    assert(ur != vr);
                    assert(!g.matrix[ur].fast_contain(vr));
                }
            }
        if (failed)
            std::cout << "oops" << std::endl;
        assert(!failed);
    }
};

void union_with_sets(
    const std::vector<int>& ind, const std::vector<bitset>& bv, bitset& result)
{
    for (auto i : ind)
        result.union_with(bv[i]);
}

bool intersect_vec_bs_p(const std::vector<int> vec, const bitset& bs)
{
    return std::any_of(
        begin(vec), end(vec), [&](int up) { return bs.fast_contain(up); });
};

// pick which partition will be our w
int pick_partition(const dense_graph& g, const bitset& clq, bitset& util_set,
    const std::vector<std::vector<int>>& partitions,
    const std::vector<int>& revmap)
{
    auto numcovered = [&](int v) {
        util_set.clear();
        union_with_sets(partitions[revmap[v]], g.origmatrix, util_set);
        // std::cout << "N(" << v << ") =" << util_set << "\n";
        return std::count_if(begin(clq), end(clq), [&](int u) {
            return u != v
                && intersect_vec_bs_p(partitions[revmap[u]], util_set);
        });
    };

    std::vector<int> nc(g.capacity());
    for (auto v : clq) {
        nc[v] = numcovered(v);
        // std::cout << "nc[" << v << "] = " << nc[v] << "\n";
    }
    return *std::min_element(begin(clq), end(clq), [&](int u, int v) {
        return std::pair<int, size_t>{-nc[u], partitions[revmap[u]].size()}
        < std::pair<int, size_t>{-nc[v], partitions[revmap[v]].size()};
    });
}

// pick a neighbor of w to remove from each partition of clq
void update_partitions(const dense_graph& g,
    std::vector<std::vector<int>>& partitions,
    std::vector<int>& covered_neighbors, const bitset& clq, int w,
    const std::vector<int>& revmap, bitset& util_set)
{
    util_set.clear();
    union_with_sets(partitions[revmap[w]], g.origmatrix, util_set);
    covered_neighbors.clear();
    for (auto u : clq) {
        auto i = std::find_if(begin(partitions[revmap[u]]),
            end(partitions[revmap[u]]),
            [&](int x) { return util_set.fast_contain(x); });
        if (i != end(partitions[revmap[u]])) {
            int x = *i;
            if (partitions[revmap[u]].size() > 1)
                partitions[revmap[u]].erase(i);
            covered_neighbors.push_back(x);
        }
    }
}

// void remap_constraint(
//     bitset& bs, bitset& util_set, const std::vector<int>& vertex_map)
// {
//     util_set.clear();
//     for (auto v : bs)
//         util_set.fast_add(vertex_map[v]);
//     bs.copy(util_set);
// }
void remap_constraint(std::vector<int>& bs, const std::vector<int>& vertex_map)
{
    for (auto vpt{begin(bs)}; vpt != end(bs); ++vpt)
        *vpt = vertex_map[*vpt];
    std::sort(begin(bs), end(bs));
    bs.erase(std::unique(begin(bs), end(bs)), end(bs));
}

cons_base* post_gc_constraint(Solver& s, dense_graph& g,
    boost::optional<fillin_info> fillin, const varmap& vars,
    const std::vector<indset_constraint>& isconses,
    const std::vector<int>& vertex_map, const options& opt, statistics& stat)
{
    // try {
    if (g.size() > 0) {
        auto cons = new gc_constraint(s, g, fillin, vars, isconses, opt, stat);
        s.addConstraint(cons);
        for (auto& c : cons->isconses) {
            // remap_constraint(c.bvs, cons->util_set, vertex_map);
            remap_constraint(c.vs, vertex_map);
        }
        return cons;
    }
    // } catch (minicsp::unsat &u) {
    return NULL;
}
} // namespace gc
