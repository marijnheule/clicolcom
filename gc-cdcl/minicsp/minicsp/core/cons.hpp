/*************************************************************************
minicsp

Copyright 2010--2014 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/

#ifndef __MINICSP_CONS_HPP
#define __MINICSP_CONS_HPP

#include "solver.hpp"
#include <vector>
#include <climits>

namespace minicsp {

/* Arithmetic relations */

// v1 == v2 + c
void post_eq(Solver& s, cspvar v1, cspvar v2, int c);

/* x == -y + c */
void post_neg(Solver &s, cspvar x, cspvar y, int c);

// v1 != v2 + c
void post_neq(Solver& s, cspvar v1, cspvar v2, int c);

// v1 <= v2 + c
void post_leq(Solver& s, cspvar v1, cspvar v2, int c);
// v1 < v2 + c
void post_less(Solver& s, cspvar v1, cspvar v2, int c);

/* v1 == v2 + c <=> b */
void post_eq_re(Solver &s, cspvar v1, cspvar v2, int c, cspvar b);
void post_eq_re(Solver &s, cspvar v1, cspvar v2, int c, Lit b);

/* v1 != v2 + c <=> b */
void post_neq_re(Solver &s, cspvar v1, cspvar v2, int c, cspvar b);
void post_neq_re(Solver &s, cspvar v1, cspvar v2, int c, Lit b);

/* v1 <= v2 + c <=> b */
void post_leq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_leq_re(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 < v2 + c <=> b */
void post_less_re(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_less_re(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 >= v2 + c <=> b */
void post_geq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_geq_re(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 > v2 + c <=> b */
void post_gt_re(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_gt_re(Solver &s, cspvar x, cspvar y, int c, Lit b);

/* right implication version of the above */
/* v1 <= v2 + c => b */
void post_leq_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_leq_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 < v2 + c => b */
void post_less_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_less_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 >= v2 + c => b */
void post_geq_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_geq_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 > v2 + c => b */
void post_gt_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_gt_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b);

/* left implication versions */
/* v1 <= v2 + c <= b */
void post_leq_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_leq_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 < v2 + c <= b */
void post_less_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_less_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 >= v2 + c <= b */
void post_geq_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_geq_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b);
/* v1 > v2 + c <= b */
void post_gt_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b);
void post_gt_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b);

/* |x| = y + c */
void post_abs(Solver& s, cspvar v1, cspvar v2, int c);

// sum coeff[i]*vars[i] + c <= 0
void post_lin_leq(Solver &s, std::vector<cspvar> const& vars,
                   std::vector<int> const &coeff, int c);
void post_lin_less(Solver &s, std::vector<cspvar> const& vars,
                    std::vector<int> const &coeff, int c);

// sum coeff[i]*vars[i] + c <= 0 implies b=1. This is NOT flatzinc's
// reified constraint!
void post_lin_leq_right_imp_re(Solver &s, std::vector<cspvar> const&vars,
                               std::vector<int> const &coeff,
                               int c, cspvar b);
void post_lin_less_right_imp_re(Solver &s, std::vector<cspvar> const&vars,
                                std::vector<int> const &coeff,
                                int c, cspvar b);

// b implies sum coeff[i]*vars[i] + c <= 0. Again, not flatzinc's
// constraint
void post_lin_leq_left_imp_re(Solver &s,
                              std::vector<cspvar> const&vars,
                              std::vector<int> const &coeff,
                              int c,
                              cspvar b);
void post_lin_less_left_imp_re(Solver &s,
                               std::vector<cspvar> const&vars,
                               std::vector<int> const &coeff,
                               int c,
                               cspvar b);

// this is flatzinc reified linear inequality
void post_lin_leq_iff_re(Solver &s, std::vector<cspvar> const& vars,
                          std::vector<int> const& coeff,
                          int c, cspvar b);
void post_lin_less_iff_re(Solver &s, std::vector<cspvar> const& vars,
                           std::vector<int> const& coeff,
                           int c, cspvar b);

/* linear equality: just <= and >=. It is np-hard to propagate anyway,
   although in the future it might be worth implementing the gac
   propagator for small coefficients */

// sum coeff[i]*vars[i] + c = 0
void post_lin_eq(Solver &s,
                 std::vector<cspvar> const& vars,
                 std::vector<int> const& coeff,
                 int c);

// sum coeff[i]*vars[i] + c = 0 implies b = 1
void post_lin_eq_right_imp_re(Solver &s,
                              std::vector<cspvar> const& vars,
                              std::vector<int> const& coeff,
                              int c,
                              cspvar b);
// b = 1 implies sum coeff[i]*vars[i] + c = 0
void post_lin_eq_left_imp_re(Solver &s,
                             std::vector<cspvar> const& vars,
                             std::vector<int> const& coeff,
                             int c,
                             cspvar b);
// b=1 iff sum coeff[i]*vars[i] + c = 0
void post_lin_eq_iff_re(Solver &s,
                        std::vector<cspvar> const& vars,
                        std::vector<int> const& coeff,
                        int c,
                        cspvar b);

/* linear inequality: L != 0 */
void post_lin_neq(Solver &s,
                  std::vector<cspvar> const& vars,
                  std::vector<int> const& coeff,
                  int c);

/* L != 0 => b = 1 */
void post_lin_neq_right_imp_re(Solver &s,
                               std::vector<cspvar> const& vars,
                               std::vector<int> const& coeff,
                               int c,
                               cspvar b);

/* b = 1 => L != 0 */
void post_lin_neq_left_imp_re(Solver &s,
                              std::vector<cspvar> const& vars,
                              std::vector<int> const& coeff,
                              int c,
                              cspvar b);

/* b = 1 <=> L != 0 */
void post_lin_neq_iff_re(Solver &s,
                         std::vector<cspvar> const& vars,
                         std::vector<int> const& coeff,
                         int c,
                         cspvar b);

/* Boolean N-ary linear inequality: Pseudo-Boolean constraints */
void post_pb(Solver& s, std::vector<Var> const& vars,
              std::vector<int> const& weights, int lb);
void post_pb(Solver& s, std::vector<cspvar> const& vars,
              std::vector<int> const& weights, int lb);

// reified versions
void post_pb_right_imp_re(Solver& s, std::vector<cspvar> const& vars,
                          std::vector<int> const& weights, int lb,
                          cspvar b);
void post_pb_left_imp_re(Solver& s, std::vector<cspvar> const& vars,
                         std::vector<int> const& weights, int lb,
                         cspvar b);
void post_pb_iff_re(Solver& s, std::vector<cspvar> const& vars,
                    std::vector<int> const& weights, int lb,
                    cspvar b);

/* Pseudo Boolean constraint with a non-Boolean var on the RHS:

   sum w[i]*v[i] + c = rhs

   where v[i] are all Boolean
 */
void post_pb(Solver& s, std::vector<Var> const& vars,
             std::vector<int> const& weights, int c, cspvar rhs);
void post_pb(Solver& s, std::vector<cspvar> const& vars,
             std::vector<int> const& weights, int c, cspvar rhs);

/* x = y*z */
void post_mult(Solver& s, cspvar x, cspvar y, cspvar z);

/* x = y/z */
void post_div(Solver& s, cspvar x, cspvar y, cspvar z);

/* x = min(y,z) */
void post_min(Solver &s, cspvar x, cspvar y, cspvar z);
/* x = max(y,z) */
void post_max(Solver &s, cspvar x, cspvar y, cspvar z);

/* Element: R = X[I-offset]

   offset is 0 by default for normal 0-based indexing, but the
   flatzinc frontend uses 1 for its own 1-based indexing */
void post_element(Solver &s, cspvar R, cspvar I,
                  std::vector<cspvar> const& X, int offset=0);

// alldiff. if gac is false it behaves like a clique of != and also
// detects gac disentailment, otherwise enforces gac like usual
void post_alldiff(Solver &s, std::vector<cspvar> const& vars, bool gac = true);

// atmostnvalue. The number of distinct values taken by the vector x
// is at most N

// BC propagator
void post_atmostnvalue(Solver &s, std::vector<cspvar> const& x, cspvar N);
// a propagator based on a heuristic for computing a max independent
// set (the MD propagator from the Bessiere et al cpaior 2005 paper)
void post_atmostnvalue_md(Solver &s, std::vector<cspvar> const& x, cspvar N);

// independent set constraint.
//
// The constraint \alpha(G) <= N, where G is a graph with $n$
// vertices, described by its adjacency matrix (see below). The
// variables x describe the upper diagonal matrix X s.t. X + X^T is
// the adjacency matrix of G.  So the first n-1 variables are the
// edges between vertex 1 and vertices 2...n, the second n-2 are the
// edges between vertex 2 and vertices 3...n and so on.
//
// Since this constraint is NP-hard even to check, it only gives a
// lower bound even for complete assignments.
void post_independent_set(Solver &s, size_t n, std::vector<Var> const& x,
                          cspvar N);

/* Regular
 *
 * the sequence x[0]...x[n] is a string accepted by the automaton (nfa
 * or dfa, does not matter) with transition function d, initial state
 * q0 and accepting states F. The alphabet S is implicitly defined by
 * the transition function.
 *
 * Since we can do NFA, we accept non-deterministic transition
 * functions but there is no way to specify e-transitions, so the nfa
 * has to be converted to e-free.
 *
 * Following the minizinc specification, state 0 is implicitly a
 * rejecting sink.
 */
namespace regular {
  struct transition {
    size_t q0;
    int s;
    size_t q1;
    transition(size_t pq0, int ps, size_t pq1) : q0(pq0), s(ps), q1(pq1) {}
  };

  struct automaton {
    std::vector<transition> d;
    int q0;
    std::set<int> F;

    automaton(std::vector<transition> const& pd, int pq0,
              std::set<int> pF) : d(pd), q0(pq0), F(pF) {}
  };
}

void post_regular(Solver& s, std::vector<cspvar> const& x,
                  regular::automaton const& aut,
                  bool gac = true);

/* Cumulative constraint: holds if a set of tasks 0..n-1 is scheduled
   so that task i starts as s[i], has duration d[i] and requires r[i]
   of a given resource, and the resource requirements of all tasks
   running at any given time is less than c. */
void post_cumulative(Solver& s, std::vector<cspvar> const& start,
                     std::vector<cspvar> const& dur,
                     std::vector<cspvar> const& req,
                     cspvar cap);


/* Table constraints. Constraints defined in extension, either as a
   set of allowed tuples (positive table) or a set of disallowed
   tuples (negative table).
 */

// exception thrown for invalid tuplesets
struct non_table{};

// tuples can contain this constant instead of an actual value. It
// means "*", which matches any value. So a positive (negative) tuple
// <1,2,*> means the constraint is satisfied (falsified) by assigning
// x1 = 1, x2 = 2, regardless of the assignment we make to x3.
static const int STAR_CONSTANT = INT_MAX;

void post_positive_table(Solver &s, std::vector<cspvar> const& x,
                         std::vector< std::vector<int> > const& tuples);
void post_negative_table(Solver &s, std::vector<cspvar> const& x,
                         std::vector< std::vector<int> > const& tuples);

/* lex ordering constraints: given two vectors x,y of equal length, it
   holds if the string x is lexicographically leq (resp, less) than
   the string y. Declaratively lex_leq(x,y) <=> x[0] <= y[0] && (x[0]
   == y[0] -> lex_leq(x[1:n], y[1:n])). lex_less(x,y) <=> lex_leq(x,y)
   && x != y. */
void post_lex_leq(Solver &s, std::vector<cspvar> const& x,
                  std::vector<cspvar> const& y);
void post_lex_less(Solver &s, std::vector<cspvar> const& x,
                   std::vector<cspvar> const& y);

} // namespace minicsp

#endif
