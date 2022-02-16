#include "rewriter.hpp"
#include "utils.hpp"

#define D                                                                      \
    if (1)                                                                     \
        ;                                                                      \
    else
#define DOUT                                                                   \
    if (1)                                                                     \
        ;                                                                      \
    else                                                                       \
        std::cout

namespace gc
{

minicsp::Solver::clause_callback_result_t rewriter::rewrite(
    vec<minicsp::Lit>& clause, int btlvl)
{
    using minicsp::domevent;

    ++numruns;
    if (s.trace)
        std::cout << "Rewriting clause #" << numruns << " = "
                  << print(s, &clause) << "\n";
    int ncol = c->ub - 1;
    partitions_eq.resize(ncol);
    partitions_neq.resize(ncol);
    for (auto& p : partitions_eq)
        p.clear();
    for (auto& p : partitions_neq)
        p.clear();
    neq_counts.clear();
    neq_sums.clear();
    neq_counts.resize(g.capacity());
    neq_sums.resize(g.capacity(), ncol * (ncol - 1) / 2);
    active_partitions.clear();
    seen.clear();
    seen.resize(s.nVars());

    for (auto l : clause)
        if (s.event(l).type != domevent::NONE)
            to_resolve.push_back(l);
    bool modified = erase_if(clause,
        [&](minicsp::Lit l) { return s.event(l).type != domevent::NONE; });

    while (!to_resolve.empty()) {
        auto l = to_resolve.back();
        to_resolve.pop_back();
        seen[var(l)] = true;
        DOUT << "Processing " << minicsp::lit_printer(s, l) << " var " << var(l)
             << std::endl;
        if (s.varLevel(l) == 0)
            continue;
        auto d = s.event(l);
        // the literals are false, so we have to file their literals
        // in the opposite sense of the event type: a false eq literal
        // describes a vertex that is *not* in the partition event.d,
        // a false neq describes a vertex that *is* in the partition
        // and so on
        switch (d.type) {
        case domevent::NONE:
            DOUT << "Using " << minicsp::lit_printer(s, l) << " as is\n";
            clause.push(l);
            break;
        case domevent::NEQ:
            DOUT << "Using " << minicsp::lit_printer(s, l) << " in partition\n";
            if (partitions_eq[d.d].empty())
                active_partitions.push_back(d.d);
            partitions_eq[d.d].push_back(d.x.id());
            break;
        case domevent::EQ:
        case domevent::LEQ:
        case domevent::GEQ: {
            auto r = s.varReason(l);
            if (!r) {
                assert(0);
            }
            for (auto o : *r)
                if (o != l && !seen[var(o)])
                    to_resolve.push_back(o);
            DOUT << minicsp::lit_printer(s, l) << " forced by " << print(s, r)
                 << "\n";
        } break;
        }
    }

    D
    {
        DOUT << "active partitions = " << print_container<std::vector<int>>{active_partitions}
             << "\n";
        for (int i = 0; i != ncol; ++i)
            DOUT << "color " << i
                 << " eq = " << print_container<std::vector<int>>{partitions_eq[i]}
                 << " neq = " << print_container<std::vector<int>>{partitions_neq[i]} << "\n";
    }

    DOUT << std::endl;

    for (int i = 0; i != ncol; ++i) {
        auto& eq = partitions_eq[i];
        auto& neq = partitions_neq[i];
        if (eq.empty())
            continue;
        int rep = eq[0];
        if (eq.size() > 1) {
            for (auto v : eq) {
                if (v == rep)
                    continue;
                assert(evars[v][rep] != minicsp::var_Undef);
                DOUT << "Adding "
                     << minicsp::lit_printer(s, ~minicsp::Lit(evars[v][rep]))
                     << "\n";
                clause.push(~minicsp::Lit(evars[v][rep]));
            }
        }
        for (auto v : neq) {
            if (evars[v][rep] != minicsp::var_Undef) {
                DOUT << "Adding "
                     << minicsp::lit_printer(s, minicsp::Lit(evars[v][rep]))
                     << "\n";
                clause.push(minicsp::Lit(evars[v][rep]));
            }
        }
    }
    for (size_t i = 0; i != active_partitions.size(); ++i)
        for (size_t j = i + 1; j != active_partitions.size(); ++j) {
            auto p0 = active_partitions[i], p1 = active_partitions[j];
            auto r0 = partitions_eq[p0][0], r1 = partitions_eq[p1][0];
            if (evars[r0][r1] != minicsp::var_Undef) {
                DOUT << "Adding "
                     << minicsp::lit_printer(s, minicsp::Lit(evars[r0][r1]))
                     << "\n";
                if (s.varLevel(evars[r0][r1]) == 0) {
                    DOUT << "\tlevel 0, ignoring\n";
                } else
                    clause.push(minicsp::Lit(evars[r0][r1]));
            }
        }

    if (modified) {
        if (s.trace)
            std::cout << "Rewrote to " << print(s, &clause) << std::endl;
        return minicsp::Solver::CCB_MODIFIED;
    } else {
        DOUT << "Not modified\n";
        return minicsp::Solver::CCB_OK;
    }
}

} // namespace gc
