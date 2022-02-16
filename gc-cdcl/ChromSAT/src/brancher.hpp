#ifndef GC_BRANCHER_HPP
#define GC_BRANCHER_HPP

#include "graph.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "utils.hpp"

#include <iomanip>

namespace gc
{

struct Brancher {
    minicsp::Solver& s;
    dense_graph& g;
    dense_graph& fg;
    const varmap& evars;
    const std::vector<minicsp::cspvar>& xvars;
    cons_base& constraint;
    const options& opt;

    int64_t numdecisions{0}, numchoices{0};

    std::mt19937 random_generator;

    Brancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg,
        const varmap& evars, const std::vector<minicsp::cspvar>& xvars,
        cons_base& constraint, const options& opt)
        : s(s)
        , g(g)
        , fg(fg)
        , evars(evars)
        , xvars(xvars)
        , constraint(constraint)
        , opt(opt)
    {
        random_generator.seed(opt.seed);
        // std::cout << random_generator() << " " << random_generator() << " "
        // << random_generator() << " " << std::endl;
    }
    virtual ~Brancher() {}

    void use()
    {
        s.varbranch = minicsp::VAR_USER;
        s.user_brancher = [this](std::vector<minicsp::Lit>& cand) {
            constraint.sync_graph();

            select_candidates(cand);
            ++numdecisions;
            numchoices += cand.size();
        };
    }

    virtual void select_candidates(std::vector<minicsp::Lit>& cand) = 0;
};

struct VSIDSBrancher : public Brancher {
    bitset util_set, low_degree;

    struct evarinfo_t {
        int u{-1}, v{-1};
    };
    std::vector<evarinfo_t> evarinfo;
    std::vector<minicsp::Var> removed;

    VSIDSBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , low_degree(0, g.capacity() - 1, bitset::empt)
    {
        for (auto u : g.nodes)
            for (auto v : g.nodes) {
                if (v <= u)
                    continue;
                if (evars[u][v] == minicsp::var_Undef)
                    continue;
                minicsp::Var var = evars[u][v];
                if (static_cast<size_t>(var) >= evarinfo.size())
                    evarinfo.resize(var + 1);
                evarinfo[var] = {u, v};
            }
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        if (!opt.branching_low_degree) {
            // plain VSIDS in this case
            auto& heap = s.vsids_heap();
            minicsp::Var next;
            do {
                next = minicsp::var_Undef;
                if (heap.empty())
                    break;
                next = heap.removeMin();
                if (s.value(next) != minicsp::l_Undef)
                    next = minicsp::var_Undef;
            } while (next == minicsp::var_Undef);
            if (next != minicsp::var_Undef)
                cand.push_back(minicsp::Lit(next));
            return;
        }

        // size_t instead of int because it gets compared to a size()
        size_t elb = std::max(*constraint.lastlb, constraint.bestlb);

        low_degree.clear();
        for (auto v : g.nodes) {
            util_set.copy(g.matrix[v]);
            util_set.intersect_with(g.nodeset);
            if (util_set.size() < elb - 1)
                low_degree.fast_add(v);
        }

        minicsp::Var next{minicsp::var_Undef};
        auto& heap = s.vsids_heap();
        removed.clear();
        do {
            next = minicsp::var_Undef;
            if (heap.empty())
                break;
            next = heap.removeMin();
            if (s.value(next) != minicsp::l_Undef) {
                next = minicsp::var_Undef;
                continue;
            }
            auto event = s.event(minicsp::Lit(next));
            if (event.type != minicsp::domevent::NONE) {
                if (low_degree.fast_contain(event.x.id())) {
                    removed.push_back(next);
                    next = minicsp::var_Undef;
                    continue;
                }
            } else {
                auto info = evarinfo[next];
                if (low_degree.fast_contain(info.u)
                    || low_degree.fast_contain(info.v)) {
                    removed.push_back(next);
                    next = minicsp::var_Undef;
                    continue;
                }
            }
        } while (next == minicsp::var_Undef);

        for (auto var : removed) {
            heap.insert(var);
        }

        if (next != minicsp::var_Undef)
            cand.push_back(minicsp::Lit(next));

        if (heap.empty()) {
            for (int i = 0; i != s.nVars(); ++i)
                assert(s.value(i) != minicsp::l_Undef);
        }
    }
};

struct VSIDSPhaseBrancher : public Brancher {
    VSIDSBrancher vsids;
    std::vector<int> candcopy;

    int D{1}, N{1};

    VSIDSPhaseBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt, int D, int N)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , vsids(s, g, fg, evars, xvars, constraint, opt)
        , D(D)
        , N(N)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        vsids.select_candidates(cand);
        assert(cand.size() <= 1);
        if (cand.empty())
            return;
        auto l = cand[0];

        auto x = var(l);
        auto event = s.event(l);

        // if VSIDS has chosen a color variable, do nothing
        if (event.type != minicsp::domevent::NONE)
            return;

        auto info = vsids.evarinfo[x];
        auto u = info.u;
        auto v = info.v;

        vsids.util_set.copy(g.matrix[u]);
        vsids.util_set.intersect_with(g.matrix[v]);
        vsids.util_set.intersect_with(g.nodeset);

        auto inter_size = vsids.util_set.size();

        vsids.util_set.copy(g.matrix[u]);
        vsids.util_set.union_with(g.matrix[v]);
        vsids.util_set.intersect_with(g.nodeset);

        auto union_size = vsids.util_set.size();

        if (inter_size * D > union_size * N)
            cand[0] = minicsp::Lit(x);
        else
            cand[0] = ~minicsp::Lit(x);
    }
};

// makes sure when choosing e_i,j that at least of i,j is in a maximal
// clique
struct VSIDSCliqueBrancher : public Brancher {
    VSIDSBrancher vsids;

    VSIDSCliqueBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , vsids(s, g, fg, evars, xvars, constraint, opt)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        // all vertices that appear in a maximal clique
        auto& cf = constraint.cf;
        auto& util_set = vsids.util_set;
        util_set.clear();
        int maxclqsz
            = *std::max_element(begin(cf.clique_sz), end(cf.clique_sz));
        for (size_t i = 0, iend = cf.clique_sz.size(); i != iend; ++i) {
            if (cf.clique_sz[i] == maxclqsz)
                util_set.union_with(cf.cliques[i]);
        }

        // choose e_ij so that at least one of i,j is in a maximal
        // clique
        auto& heap = s.vsids_heap();
        auto& removed = vsids.removed;
        removed.clear();
        minicsp::Var next;
        do {
            next = minicsp::var_Undef;
            if (heap.empty())
                break;
            next = heap.removeMin();
            if (s.value(next) != minicsp::l_Undef) {
                next = minicsp::var_Undef;
                continue;
            }
            if (s.event(minicsp::Lit(next)).type != minicsp::domevent::NONE)
                break;
            auto varinfo = vsids.evarinfo[next];
            if (!util_set.fast_contain(varinfo.u)
                && !util_set.fast_contain(varinfo.v)) {
                removed.push_back(next);
                next = minicsp::var_Undef;
            }
        } while (next == minicsp::var_Undef);
        if (next != minicsp::var_Undef)
            cand.push_back(minicsp::Lit(next));
        for (auto var : removed)
            heap.insert(var);
        return;
    }
};

struct VSIDSColorBrancher : public Brancher {
    VSIDSColorBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
    {
        for (auto& uv : evars.vars)
            for (auto ep : uv.second) {
                auto e = ep.second;
                if (e != minicsp::var_Undef)
                    s.setDecisionVar(e, false);
            }
        for (auto x : xvars)
            for (int i = x.omin(s); i <= x.omax(s); ++i) {
                auto v = x.leqi(s, i);
                if (v != minicsp::var_Undef)
                    continue;
                s.setDecisionVar(v, false);
            }
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        auto& heap = s.vsids_heap();
        minicsp::Var next;
        do {
            next = minicsp::var_Undef;
            if (heap.empty())
                break;
            next = heap.removeMin();
            if (s.value(next) != minicsp::l_Undef)
                next = minicsp::var_Undef;
            auto event = s.event(minicsp::Lit(next));
            if (event.type != minicsp::domevent::EQ) {
                next = minicsp::var_Undef;
                continue;
            }
        } while (next == minicsp::var_Undef);
        if (next != minicsp::var_Undef) {
            auto l = minicsp::Lit(next);
            cand.push_back(l);
        }
        return;
    }
};

struct BrelazBrancher : public Brancher {
    std::vector<int> mindom;
    // a maximal clique in 3 forms: vector, bitset and remapped to the
    // representatives
    std::vector<int> clique;
    bitset clique_bs;

    std::vector<int> clqorder;

    // temp
    bitset util_set, util_set2;

    // vertices that are ignored because they have low degree
    std::vector<int> low_degree;

    int nlit{0};
    std::vector<minicsp::Lit> level_0_choices;

    BrelazBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , clique_bs(0, g.capacity() - 1, bitset::empt)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , util_set2(0, g.capacity() - 1, bitset::empt)
    {
        auto& cf = constraint.cf;
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        auto& clq = cf.cliques[maxidx];
        for (auto v : clq)
            clique.push_back(v);
    }

    void select_candidates_xvars(std::vector<minicsp::Lit>& cand)
    {
        for (auto v : clique) {
            auto x = xvars[v];
            if (x.domsize(s) != 1) {
                cand.clear();
                cand.push_back(x.e_eq(s, x.min(s)));
                return;
            }
        }

        // size_t instead of int because it gets compared to a size()
        size_t elb = std::max(*constraint.lastlb, constraint.bestlb);

        int mind{-1};
        mindom.clear();
        low_degree.clear();
        for (auto v : g.nodes) {
            auto x = xvars[v];
            auto xd = x.domsize(s);
            if (xd == 1)
                continue;
            if (opt.branching_low_degree) {
                util_set.copy(g.matrix[v]);
                util_set.intersect_with(g.nodeset);
                if (util_set.size() < elb) {
                    low_degree.push_back(v);
                    continue;
                }
            }
            if (mindom.empty() || xd < mind) {
                mindom.clear();
                mind = x.domsize(s);
                mindom.push_back(v);
            } else if (xd == mind)
                mindom.push_back(v);
        }

        if (mindom.empty()) {
            if (!low_degree.empty()) {
                auto x = xvars[low_degree.back()];
                cand.clear();
                cand.push_back(x.e_eq(s, x.min(s)));
                return;
            } else
                return;
        }

        // int tiedv = mindom.size();
        int maxdv{-1};
        int maxd{-1};
        if (mindom.size() > 1) {
            for (auto v : mindom) {
                util_set.copy(g.matrix[v]);
                util_set.intersect_with(g.nodeset);
                int deg = util_set.size();
                if (deg > maxd) {
                    maxd = deg;
                    maxdv = v;
                }
            }
        } else
            maxdv = mindom[0];
        auto x = xvars[maxdv];
        // std::cout << "Branching on vertex " << maxdv
        //           << " domsize = " << x.domsize(s) << " degree = " << maxd
        //           << " domsize ties = " << mindom.size() << "\n";
        cand.clear();
        cand.push_back(x.e_eq(s, x.min(s)));
    }

    void select_candidates_evars(std::vector<minicsp::Lit>& cand)
    {
        auto& cf = constraint.cf;
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        auto& clq = cf.cliques[maxidx];

        clique_bs.copy(clq);
        if (clique_bs.size() == g.nodeset.size())
            return;

        // std::cout << clq.size() << " " << clique_bs;

        util_set.copy(clique_bs);
        util_set.intersect_with(g.nodeset);
        int clqsize = util_set.size();

        assert(util_set.size() == clique_bs.size());

        int maxc{-1}, maxd{-1}, maxv{-1};
        for (auto v : g.nodes) {
            if (clique_bs.fast_contain(v))
                continue;

            // number of neighboring colors == intersection of
            // neighborhood with clique
            util_set.copy(g.matrix[v]);
            util_set.intersect_with(clique_bs);
            int vc = util_set.size();

            assert(vc < clqsize);

            if (vc < maxc || vc == clqsize)
                continue;

            // degree == neighborhood setminus clique (restricted to
            // current graph)
            util_set.copy(g.matrix[v]);
            util_set.setminus_with(clique_bs);
            util_set.intersect_with(g.nodeset);
            int vd = util_set.size();

            // max neighboring colors, tie breaking by degree
            if ((vc > maxc && vc < clqsize) || (vc == maxc && vd > maxd)) {
                maxv = v;
                maxc = vc;
                maxd = vd;
            }
        }

        assert(maxv >= 0);

        util_set.copy(clique_bs);
        util_set.setminus_with(g.matrix[maxv]);
        int u = (random_generator() % 2 ? util_set.min() : util_set.max());

        // std::cout << util_set.size() << std::endl;
        minicsp::Var evar{minicsp::var_Undef};
        if (opt.fillin) {
            util_set.clear();
            for (auto up : g.partition[u])
                util_set.fast_add(up);
            // find an actual var that can merge u and maxv
            for (auto vp : g.partition[maxv]) {
                util_set2.copy(util_set);
                util_set2.intersect_with(constraint.fg.matrix[vp]);
                if (!util_set2.empty()) {
                    evar = evars[vp][util_set2.min()];
                    break;
                }
            }
            assert(evar != minicsp::var_Undef);
        } else {
            evar = evars[u][maxv];

            // std::cout << " - " << maxc << " - (" << u << "=" << maxv << ")"
            //           << std::endl;
        }
        assert(evar != minicsp::var_Undef);
        cand.push_back(minicsp::Lit(evar));
    }

    void select_candidates_cvars(std::vector<minicsp::Lit>& cand)
    {

        if (s.decisionLevel() <= 1 and nlit == s.nVars()) {

            // std::cout << "Use previous branching decisions\n";

            for (auto l : level_0_choices) {
                cand.push_back(l);
            }

            return;
        }

        auto& cf = constraint.cf;

        // if (s.decisionLevel() <= 1 and nlit != s.nVars()) {
        //     std::cout << "Number of unit literal changed! (" << s.nVars()
        //               << " vars)\n";
        // }

        clqorder.clear();
        for (int i = 0; i < cf.num_cliques; ++i)
            clqorder.push_back(i);
        std::sort(begin(clqorder), end(clqorder),
            [&](int a, int b) { return cf.clique_sz[a] > cf.clique_sz[b]; });

        // std::cout << "level = " << s.decisionLevel() << ", " << g.size() << " nodes, " << g.count_edges()
        //           << " edges, max clique = " << cf.clique_sz[clqorder[0]]
        //           << std::endl;

        size_t maxidx{0};
        // explore the cliques by decreasing size until we find one with a
        // missing chord
        while (maxidx < clqorder.size()) {

            clique_bs.copy(cf.cliques[clqorder[maxidx]]);
            // int sz_before = clique_bs.size();
            // clique_bs.intersect_with(g.nodeset);
            // assert(sz_before == clique_bs.size());

            // std::cout << "try " << clique_bs << std::endl;

            bool missing_chord = false;

            // try to find a vertex v of the
            // clique s.t. (Chord(v) & V) \ N(v)
            // is non-empty
            for (auto v : clique_bs) {
                util_set.copy(fg.origmatrix[v]);
                util_set.intersect_with(g.nodeset);

                if (!g.matrix[v].includes(util_set)) {

                    gc::bitset mchords(0, g.capacity() - 1, gc::bitset::empt);
                    mchords.copy(util_set);
                    mchords.setminus_with(g.matrix[v]);
                    // std::cout << " -> missing "
                    //              "chords! "
                    //           << mchords << std::endl;

                    missing_chord = true; // ok we have found one
                    break;
                }
            }
            if (missing_chord)
                break;
            else {
                // move to the next clique
                ++maxidx;
            }
        }

        if (maxidx >= clqorder.size())
            return;

        clique_bs.copy(cf.cliques[clqorder[maxidx]]);
        // clique_bs.intersect_with(g.nodeset);
        util_set2.clear(); // we are going to store all the potential chords in
        // util_set2

        for (auto v : clique_bs) {
            util_set.copy(fg.origmatrix[v]);
            util_set.intersect_with(g.nodeset);
            util_set.setminus_with(g.matrix[v]);
            util_set2.union_with(util_set);
        }

        // now pick the node with highest saturation (w.r.t. this clique)
        int maxc{-1}, maxd{-1}, maxv{-1};
        for (auto v : util_set2) {
            assert(!clique_bs.fast_contain(v));

            // number of neighboring colors == intersection of
            // neighborhood with clique
            util_set.copy(g.matrix[v]);
            util_set.intersect_with(clique_bs);
            int vc = util_set.size();
            if (vc < maxc)
                continue;

            // degree == neighborhood setminus clique (restricted to
            // current graph)
            util_set.copy(g.matrix[v]);
            util_set.setminus_with(clique_bs);
            util_set.intersect_with(g.nodeset);
            int vd = util_set.size();

            // max neighboring colors, tie breaking by degree
            if (vc > maxc || (vc == maxc && vd > maxd)) {
                maxv = v;
                maxc = vc;
                maxd = vd;
            }
        }

        assert(maxv >= 0);

        util_set.copy(clique_bs);
        util_set.intersect_with(fg.origmatrix[maxv]);
        util_set.setminus_with(g.matrix[maxv]);

        // int u = util_set.min();

        for (auto u : util_set) {
            minicsp::Var evar{minicsp::var_Undef};
            evar = evars[u][maxv];
            assert(evar != minicsp::var_Undef);
            cand.push_back(minicsp::Lit(evar));
        }

        if (s.decisionLevel() <= 1 and nlit != s.nVars()) {
            nlit = s.nVars();
            level_0_choices.clear();
            for (auto l : cand)
                level_0_choices.push_back(l);
        }
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        // if (s.decisionLevel() <= 1 and nlit == s.nVars()) {
        //
        //     // std::cout << "Use previous branching decisions\n";
        //
        //     for (auto l : level_0_choices) {
        //         cand.push_back(l);
        //     }
        //
        //     return;
        // }

        if (opt.xvars)
            select_candidates_xvars(cand);
        else if (opt.fillin)
            select_candidates_cvars(cand);
        else
            select_candidates_evars(cand);

        // if (!opt.xvars and s.decisionLevel() <= 1 and nlit != s.nVars()) {
        //     nlit = s.nVars();
        //     level_0_choices.clear();
        //     for (auto l : cand)
        //         level_0_choices.push_back(l);
        // }
    }
};

template <int N, int D> struct EdgeBrancher : public Brancher {
    using Brancher::Brancher;
    std::vector<edge> e_cand;
    std::vector<int> nodes;
    bitset neighbors_u;
    bitset neighbors_v;
    bitset counter;
    int max_tied;

    EdgeBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , neighbors_u(0, g.capacity() - 1, bitset::empt)
        , neighbors_v(0, g.capacity() - 1, bitset::empt)
        , counter(0, g.capacity() - 1, bitset::empt)
        , max_tied(g.capacity() * g.capacity())
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand) = 0;

    void select_branch(std::vector<minicsp::Lit>& cand)
    {
        for (auto e : e_cand) {
            if (N < 0)
                cand.push_back(minicsp::Lit(evars[e.first][e.second]));
            else {
                auto u{e.first};
                auto v{e.second};

                counter.copy(g.matrix[u]);
                counter.intersect_with(g.matrix[v]);
                counter.intersect_with(g.nodeset);

                auto inter_size = counter.size();

                counter.copy(g.matrix[u]);
                counter.union_with(g.matrix[v]);
                counter.intersect_with(g.nodeset);

                auto union_size = counter.size();

                if (inter_size * D > union_size * N) {
                    cand.push_back(minicsp::Lit(evars[u][v]));
                } else {
                    cand.push_back(~minicsp::Lit(evars[u][v]));
                }
            }
        }
    }
};

template <int N, int D, typename Op>
struct PartitionBrancher : public EdgeBrancher<N, D> {
    Op op;

    PartitionBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt, Op op)
        : EdgeBrancher<N, D>(s, g, fg, evars, xvars, constraint, opt)
        , op(op)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            auto u_part_size{this->g.partition[u].size()};
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;

                auto criterion{op(u_part_size, this->g.partition[v].size())};
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

template <int N, int D, typename Op>
std::unique_ptr<PartitionBrancher<N, D, Op>> make_partition_brancher(
    minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
    const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
    const options& opt, Op op)
{
    return std::make_unique<PartitionBrancher<N, D, Op>>(
        s, g, fg, evars, xvars, constraint, opt, op);
}

template <int N, int D, typename Op>
struct DegreeBrancher : public EdgeBrancher<N, D> {
    Op op;

    DegreeBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt, Op op)
        : EdgeBrancher<N, D>(s, g, fg, evars, xvars, constraint, opt)
        , op(op)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            this->neighbors_u.copy(this->g.matrix[u]);
            this->neighbors_u.intersect_with(this->g.nodeset);
            auto u_degree{this->neighbors_u.size()};
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;

                this->neighbors_v.copy(this->g.matrix[v]);
                this->neighbors_v.intersect_with(this->g.nodeset);

                auto criterion{op(this->neighbors_v.size(), u_degree)};
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

template <int N, int D, typename Op>
std::unique_ptr<DegreeBrancher<N, D, Op>> make_degree_brancher(
    minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
    const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
    const options& opt, Op op)
{
    return std::make_unique<DegreeBrancher<N, D, Op>>(
        s, g, fg, evars, xvars, constraint, opt, op);
}

template <int N, int D> struct DegreeUnionBrancher : public EdgeBrancher<N, D> {
    // using EdgeBrancher::EdgeBrancher;

    DegreeUnionBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg, const varmap& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : EdgeBrancher<N, D>(s, g, fg, evars, xvars, constraint, opt)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            this->neighbors_u.copy(this->g.matrix[u]);
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;
                this->counter.copy(this->g.matrix[v]);
                this->counter.intersect_with(this->g.nodeset);
                int criterion = this->counter.size();
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

// computes an activity for each vertex which is increased every time
// a vertex participates in conflict resolution and decays
// exponentially. On branching, select the largest clique and merge
// the highest activity vertex outside the clique with the highest
// activity vertex in the clique
struct VertexActivityBrancher : public Brancher {
    double vtx_inc{1};
    double vtx_decay{1 / 0.95};
    std::vector<double> activity;
    bitset seen;
    bitset util_set;
    BrelazBrancher brelaz;

    VertexActivityBrancher(minicsp::Solver& s, dense_graph& g, dense_graph& fg,
        const varmap& evars, const std::vector<minicsp::cspvar>& xvars,
        cons_base& constraint, const options& opt)
        : Brancher(s, g, fg, evars, xvars, constraint, opt)
        , activity(g.capacity())
        , seen(0, g.capacity() - 1, bitset::empt)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , brelaz(s, g, fg, evars, xvars, constraint, opt)
    {
        s.use_clause_callback([this](auto& cls, int) { return clscb(cls); });
    }

    void vtxBumpActivity(int v)
    {
        activity[v] += vtx_inc;
        if (activity[v] > 1e100) {
            // Rescale:
            for (int i = 0; i != g.capacity(); i++)
                activity[i] *= 1e-100;
            vtx_inc *= 1e-100;
        }
    }

    void vtxDecayActivity() { vtx_inc *= vtx_decay; };

    minicsp::Solver::clause_callback_result_t clscb(vec<minicsp::Lit>& cls)
    {
        seen.clear();
        for (auto l : cls) {
            auto info = constraint.varinfo[var(l)];
            if (!seen.fast_contain(info.u)) {
                vtxBumpActivity(info.u);
                seen.fast_add(info.u);
            }
            if (!seen.fast_contain(info.v)) {
                vtxBumpActivity(info.v);
                seen.fast_add(info.v);
            }
        }
        vtxDecayActivity();
        return minicsp::Solver::CCB_OK;
    }

    void select_candidates(std::vector<minicsp::Lit>& cand) override
    {
        if (opt.brelaz_first && s.conflicts < 1e5)
            return brelaz.select_candidates_evars(cand);

        auto& cf = constraint.cf;
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        auto& clq = cf.cliques[maxidx];

        if (clq.size() == g.nodeset.size())
            return;

        int vtx{-1};
        double vtxact{0.01};
        size_t vtxdom{g.nodes.size()};
        for (auto v : g.nodes) {
            if (clq.fast_contain(v))
                continue;
            double vact = std::max(activity[v], 0.01);
            double vdom = [&]() {
                if (opt.branching == options::VERTEX_ACTIVITY) {
                    // skips bitset stuff if unnecessary
                    return 0u;
                } else {
                    util_set.copy(g.matrix[v]);
                    util_set.intersect_with(clq);
                    return util_set.size();
                }
                assert(0);
            }();
            switch (opt.branching) {
            case options::VERTEX_ACTIVITY:
                if (vtx < 0 || vact > vtxact) {
                    vtx = v;
                    vtxact = vact;
                }
                break;
            case options::VERTEX_DOM_OVER_ACT:
                if (vtx < 0 || vdom / vact > vtxact / vtxdom) {
                    vtx = v;
                    vtxact = vact;
                    vdom = vtxdom;
                }
                break;
            case options::VERTEX_DOM_THEN_ACT:
                if (vtx < 0 || vdom < vtxdom
                    || (vdom == vtxdom && vact > vtxact)) {
                    vtx = v;
                    vtxdom = vdom;
                    vtx = vact;
                }
                break;
            default:
                assert(0);
            }
        }
        assert(vtx >= 0);

        int utx{-1};
        for (auto u : clq) {
            if (g.matrix[vtx].fast_contain(u))
                continue;
            if (utx < 0 || activity[u] > activity[utx])
                utx = u;
        }
        assert(utx >= 0);
        auto evar = evars[utx][vtx];
        assert(evar != minicsp::var_Undef);
        if (opt.phase_saving)
            cand.push_back(minicsp::Lit(
                evar, (s.currentVarPhase(evar) == minicsp::l_False)));
        else
            cand.push_back(minicsp::Lit(evar));
    }
};

} // namespace gc

#endif
