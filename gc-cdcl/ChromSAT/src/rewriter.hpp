#ifndef __GC_REWRITER_HPP
#define __GC_REWRITER_HPP

#include "graph.hpp"
#include "minicsp/core/solver.hpp"
#include "prop.hpp"

namespace gc
{
struct rewriter {
    minicsp::Solver& s;
    const dense_graph& g;
    cons_base* c;
    const varmap& evars;
    const std::vector<minicsp::cspvar>& xvars;

    // which partitions have at least one vertex in them
    std::vector<int> active_partitions;
    // per color: array of variables in that partition, and
    // excluded from that partition
    std::vector<std::vector<int>> partitions_eq;
    std::vector<std::vector<int>> partitions_neq;
    // A clause can describe x in partition d by saying it is not
    // in any other partition. This keeps the number of partitions
    // each vertex is not in, so if it is ub-1, we can put it in
    // the remaining partition.
    std::vector<int> neq_counts;
    // this keeps sum_d - sum{d | i in partitions_neq[d]}. If
    // neq_counts[i]==ub-1 then neq_sums gives the partition it is
    // in.
    std::vector<int> neq_sums;

    std::vector<minicsp::Lit> to_resolve;
    std::vector<uint8_t> seen;

    int numruns{0};

    rewriter(minicsp::Solver& s, const dense_graph& g, cons_base* c,
        const varmap& evars, const std::vector<minicsp::cspvar>& xvars)
        : s(s)
        , g(g)
        , c(c)
        , evars(evars)
        , xvars(xvars)
    {
    }

    minicsp::Solver::clause_callback_result_t rewrite(
        vec<minicsp::Lit>& clause, int btlvl);
};
} // namespace gc

#endif
