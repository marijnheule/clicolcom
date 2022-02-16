/*****************************************************************************
minicsp

Copyright 2009--2011 George Katsirelos

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

Parts of this file were distributed under the following license as
part of MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include "solver.hpp"
#include "cons.hpp"
#include "minicsp/mtl/Sort.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

namespace minicsp {

//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters: (formerly in 'SearchParams')
    trace(false)
  , debugclauses(false)
  , learning(true)
  , restarting(true)
  , var_decay(1 / 0.95), clause_decay(1 / 0.999), random_var_freq(0.02)
  , restart_first(32), restart_inc(1.5), learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // More parameters:
    //
  , expensive_ccmin  (true)
  , polarity_mode    (polarity_false)
  , verbosity        (0)
  , phase_saving     (true)
  , solution_phase_saving(false)
  , allow_clause_dbg (true)

  , conflict_lim     (-1)

    // branching heuristics
  , varbranch(VAR_VSIDS)
  , valbranch(VAL_LEX)
    // Statistics: (formerly in 'SolverStats')
    //
  , starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
  , clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , ok               (true)
  , cla_inc          (1)
  , var_inc          (1)
  , qhead            (0)
  , simpDB_assigns   (-1)
  , simpDB_props     (0)
  , order_heap       (VarOrderLt(activity))
  , random_seed      (91648253)
  , progress_estimate(0)
  , remove_satisfied (true)

    , backtrackable_size(0)
    , backtrackable_cap(0)
    , backtrackable_space(0)
    , current_space(0L)

    , active_constraint(0L)
{
  consqs.growTo(MAX_PRIORITY+2);
  reset_queue();
  prop_queue = &consqs[0];
}


Solver::~Solver()
{
  for (int i = 0; i < learnts.size(); i++) free(learnts[i]);
  for (int i = 0; i < clauses.size(); i++) free(clauses[i]);
  for (int i = 0; i < conses.size(); ++i) conses[i]->dispose();
  for (int i = 0; i != inactive.size(); ++i) free(inactive[i]);
  for (int i = 0; i != cspvars.size(); ++i) {
    cspvar_fixed & xf = cspvars[i];
    if( xf.omax == xf.omin ) continue;
    for(int j = 0; j != (xf.omax-xf.omin+1); ++j) {
      free(xf.ps1[j]);
      free(xf.ps2[j]);
      if( j > 0 ) free(xf.ps3[j]);
      free(xf.ps4[j]);
    }
  }
  for(int i = 0; i != backtrackable_space.size(); ++i)
    free(backtrackable_space[i]);
  free(current_space);
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches   .push();          // (list for positive literal)
    watches   .push();          // (list for negative literal)
    binwatches   .push();          // (list for positive literal)
    binwatches   .push();          // (list for negative literal)

    wakes_on_lit.push();
    sched_on_lit.push();

    reason    .push({});
    assigns   .push(toInt(l_Undef));
    level     .push(-1);
    activity  .push(0);
    seen      .push(0);

    polarity    .push((char)sign);
    decision_var.push((char)dvar);

    phase.push(l_Undef);

    domevent none;
    events.push(none);
    events.push(none);
    setevent setnone;
    setevents.push(setnone);
    setevents.push(setnone);

    varnames.push_back(std::string());

    insertVarOrder(v);
    return v;
}

cspvar Solver::newCSPVar(int min, int max)
{
  assert(max - min >= 0 );

  bool unary = false;
  if( max == min ) unary = true;

  cspvar x(cspvars.size());

  cspvars.push();
  cspvarnames.push_back(std::string());

  cspvar_fixed & xf = cspvars.last();

  xf.omin = min;
  xf.omax = max;
  xf.min = min;
  xf.max = max;
  xf.dsize = max-min+1;

  reduce_var_seen.push_back(false);
  reduce_var_min.push_back(xf.omin);
  reduce_var_max.push_back(xf.omax);
  reduce_var_asgn.push_back(false);

  // the propositional encoding of the domain
  xf.firstbool = newVar();
  for(int i = 1; i != 2*xf.dsize; ++i)
    newVar();

  for(int i = xf.omin; i <= xf.omax; ++i) {
    domevent eq(x, domevent::EQ, i),
      neq(x, domevent::NEQ, i),
      leq(x, domevent::LEQ, i),
      geq(x, domevent::GEQ, i+1);
    events[ toInt( Lit(xf.eqi(i) ) ) ] = eq;
    events[ toInt( ~Lit(xf.eqi(i) ) ) ] = neq;
    events[ toInt( Lit(xf.leqi(i) ) ) ] = leq;
    events[ toInt( ~Lit(xf.leqi(i) ) ) ] = geq;
  }

  if( !unary ) {
    xf.ps1.growTo(xf.dsize);
    xf.ps2.growTo(xf.dsize);
    xf.ps3.growTo(xf.dsize);
    xf.ps4.growTo(xf.dsize);
  }

  // (x <= i) => (x <= i+1)
  // (x = i) <=> (x <= i) /\ -(x <= i-1)
  for(int i = 0; i != xf.dsize; ++i) {
    if( unary ) continue;
    vec<Lit> ps1, ps2, ps3, ps4;
    pushifdef(ps1,  ~Lit(xf.leqi(i+xf.omin)) );
    pushifdef(ps1, Lit(xf.leqi(i + 1 + xf.omin)));
    Clause *c1 = Clause_new(ps1);
    xf.ps1[i] = c1;

    ps2.push( ~Lit(xf.eqi(i+xf.omin)) );
    ps2.push( Lit(xf.leqi(i+xf.omin)) );
    Clause *c2 = Clause_new(ps2);
    xf.ps2[i] = c2;

    if( i > 0 ) {
      ps3.push( ~Lit(xf.eqi(i+xf.omin)) );
      ps3.push( ~Lit(xf.leqi(i-1+xf.omin) ) );
      Clause *c3 = Clause_new(ps3);
      xf.ps3[i] = c3;

      ps4.push( ~Lit(xf.leqi(i+xf.omin)) );
      ps4.push( Lit(xf.leqi(i-1+xf.omin)) );
      ps4.push( Lit(xf.eqi(i+xf.omin)) );
      Clause *c4 = Clause_new(ps4);
      xf.ps4[i] = c4;
    } else {
      xf.ps3[i] = INVALID_CLAUSE;
      ps4.push( ~Lit(xf.leqi(xf.omin)) );
      ps4.push( Lit(xf.eqi(xf.omin)) );
      Clause *c4 =Clause_new(ps4);
      xf.ps4[i] = c4;
    }
  }

  /* x <= omax, so this is true always. We could just simplify the
     rest of the clauses but this seems easier
   */
  if( !unary )
    x.setmax(*this, xf.omax, (Clause*)0L);

  /* Unary vars are also hacky. Immediately set eqi(min) = l_True
   */
  if( unary ) {
    uncheckedEnqueue_np(Lit(xf.firstbool+1), NO_REASON);
    uncheckedEnqueue_np(Lit(xf.firstbool), NO_REASON);
  }

  return x;
}

std::vector<cspvar> Solver::newCSPVarArray(int n, int min, int max)
{
  std::vector<cspvar> rv;
  for(int i = 0; i != n; ++i)
    rv.push_back(newCSPVar(min, max));
  return rv;
}

setvar Solver::newSetVar(int min, int max)
{
  assert(max - min >= 0 );

  setvar x(setvars.size());

  setvars.push();
  setvarnames.push_back(std::string());

  setvar_data & xd = setvars.last();

  xd.min = min;
  xd.max = max;
  xd._card = newCSPVar(0, max-min+1);

  // the propositional encoding of the domain
  xd.firstbool = newVar();
  for(int i = min+1; i != max+1; ++i)
    newVar();

  for(int i = xd.min; i <= xd.max; ++i) {
    setevent in(x, setevent::IN, i),
      ex(x, setevent::EX, i);
    setevents[ toInt( Lit(xd.ini(i) ) ) ] = in;
    setevents[ toInt( ~Lit(xd.ini(i) ) ) ] = ex;
  }

  /* post constraint sum_i ini(i) = _card */
  vector<Var> v(max-min+1);
  vector<int> w(max-min+1, 1);
  for(int i = min; i != max+1; ++i)
    v[i-min] = xd.firstbool + i - min;
  post_pb(*this, v, w, 0, xd._card);

  return x;
}

std::vector<setvar> Solver::newSetVarArray(int n, int min, int max)
{
  std::vector<setvar> rv;
  for(int i = 0; i != n; ++i)
    rv.push_back(newSetVar(min, max));
  return rv;
}

void Solver::addClause(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);

    if (!ok)
      throw unsat();
    else{
        // Check if clause is satisfied and remove false/duplicate literals:
        sort(ps);
        Lit p; int i, j;
        for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
            assert(var(ps[i]) != var_Undef );
            if (value(ps[i]) == l_True || ps[i] == ~p)
                return;
            else if (value(ps[i]) != l_False && ps[i] != p)
                ps[j++] = p = ps[i];
        }
        ps.shrink(i - j);
    }

    if (ps.size() == 0) {
      ok = false;
      throw unsat();
    } else if (ps.size() == 1){
        assert(value(ps[0]) == l_Undef);
        uncheckedEnqueue(ps[0]);
        ok = (propagate() == NULL);
        if( !ok ) throw unsat();
        return;
    }else{
        Clause* c = Clause_new(ps, false);
        clauses.push(c);
        attachClause(*c);
        if( trace )
          cout << "Added clause " << print(*this, c) << "\n";
    }
}

void Solver::addInactiveClause(Clause* c)
{
  inactive.push(c);
}

void Solver::attachClause(Clause& c) {
    assert(c.size() > 1);
    if (c.size() == 2) {
        binwatches[toInt(~c[0])].push({&c, c[1]});
        binwatches[toInt(~c[1])].push({&c, c[0]});
    } else {
        watches[toInt(~c[0])].push({&c, c[1]});
        watches[toInt(~c[1])].push({&c, c[0]});
    }
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(Clause& c) {
    assert(c.size() > 1);
    if (c.size() == 2) {
      assert(find(binwatches[toInt(~c[0])], watch{&c, lit_Undef}));
      assert(find(binwatches[toInt(~c[1])], watch{&c, lit_Undef}));
      remove(binwatches[toInt(~c[0])], watch{&c, lit_Undef});
      remove(binwatches[toInt(~c[1])], watch{&c, lit_Undef});
    } else {
      assert(find(watches[toInt(~c[0])], watch{&c, lit_Undef}));
      assert(find(watches[toInt(~c[1])], watch{&c, lit_Undef}));
      remove(watches[toInt(~c[0])], watch{&c, lit_Undef});
      remove(watches[toInt(~c[1])], watch{&c, lit_Undef});
    }
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(Clause& c) {
    detachClause(c);
    free(&c); }


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }

bool Solver::addConstraint(cons *c)
{
  conses.push(c);
  return true;
}

void Solver::wake_on_lit(Var v, cons *c, void *advice)
{
  wakes_on_lit[v].push( make_pair(c, advice) );
}

void Solver::wake_on_dom(cspvar x, cons *c, void *advice)
{
  cspvars[x._id].wake_on_dom.push( make_pair(c, advice) );
}

void Solver::wake_on_lb(cspvar x, cons *c, void *advice)
{
  cspvars[x._id].wake_on_lb.push( make_pair(c, advice) );
}

void Solver::wake_on_ub(cspvar x, cons *c, void *advice)
{
  cspvars[x._id].wake_on_ub.push( make_pair(c, advice) );
}

void Solver::wake_on_fix(cspvar x, cons *c, void *advice)
{
  cspvars[x._id].wake_on_fix.push( make_pair(c, advice) );
}

void Solver::ensure_can_schedule(cons *c)
{
  assert( c->priority >= 0 && c->priority <= MAX_PRIORITY );
  if( c->cqidx < 0 ) {
    c->cqidx = consqs.size();
    consqs.push( consqueue(c, c->priority) );
  }
}

void Solver::schedule_on_lit(Var x, cons *c)
{
  ensure_can_schedule(c);
  sched_on_lit[x].push(c->cqidx);
}

void Solver::schedule_on_dom(cspvar x, cons *c)
{
  ensure_can_schedule(c);
  cspvars[x._id].schedule_on_dom.push(c->cqidx);
}

void Solver::schedule_on_lb(cspvar x, cons *c)
{
  ensure_can_schedule(c);
  cspvars[x._id].schedule_on_lb.push(c->cqidx);
}

void Solver::schedule_on_ub(cspvar x, cons *c)
{
  ensure_can_schedule(c);
  cspvars[x._id].schedule_on_ub.push(c->cqidx);
}

void Solver::schedule_on_fix(cspvar x, cons *c)
{
  ensure_can_schedule(c);
  cspvars[x._id].schedule_on_fix.push(c->cqidx);
}

void Solver::setVarName(Var v, std::string const& name)
{
  varnames[v] = name;
}

void Solver::setCSPVarName(cspvar v, std::string const& name)
{
  cspvarnames[v._id] = name;
}

void Solver::setSetVarName(setvar v, std::string const& name)
{
  setvarnames[v._id] = name;
}

std::string const& Solver::getVarName(Var v)
{
  return varnames[v];
}

std::string const& Solver::getCSPVarName(cspvar v)
{
  return cspvarnames[v._id];
}

std::string const& Solver::getSetVarName(setvar v)
{
  return setvarnames[v._id];
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var     x  = var(trail[c]);
            if( phase_saving )
              phase[x]   = toLbool(assigns[x]);
            assigns[x] = toInt(l_Undef);
            reason[x].reset();
            insertVarOrder(x);
            domevent const &pevent = events[toInt(trail[c])];
            if( noevent(pevent) ) continue;
            cspvar_fixed& xf = cspvars[pevent.x._id];
            switch( pevent.type ) {
            case domevent::NEQ: ++xf.dsize; break;
            case domevent::LEQ: xf.max = max(xf.max, pevent.d+1); break;
            case domevent::GEQ: xf.min = min(xf.min, pevent.d-1); break;
            default: break;
            }
        }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
        memcpy(current_space, backtrackable_space[level],
               backtrackable_size);
    }
}

btptr Solver::alloc_backtrackable(unsigned size)
{
  // it is not strictly necessary that we are at level 0, but it is
  // easier to assert this until we find a use case that requires
  // otherwise
  assert(decisionLevel() == 0);
  assert(size > 0);
  btptr p;
  // make sure the allocation is sizeof(int)-aligned
  p.offset = unsigned(ceil(double(backtrackable_size)/sizeof(int)))*sizeof(int);
  backtrackable_size = p.offset+size;
  if( backtrackable_size > backtrackable_cap ) {
    backtrackable_cap = std::max(backtrackable_cap*2,
                                 backtrackable_size);
    current_space = realloc(current_space, backtrackable_cap);
    if( !current_space ) throw std::bad_alloc();

    // free cached allocations
    for(int i = 0; i != backtrackable_space.size(); ++i) {
      free(backtrackable_space[i]);
      backtrackable_space[i] = 0L;
    }
  }
  return p;
}

//=================================================================================================
// Major methods:


Lit Solver::pickBranchLitVSIDS(int polarity_mode, double random_var_freq)
{
    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (toLbool(assigns[next]) == l_Undef && decision_var[next])
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || toLbool(assigns[next]) != l_Undef || !decision_var[next])
        if (order_heap.empty()){
            next = var_Undef;
            break;
        }else
            next = order_heap.removeMin();

    if( next == var_Undef )
        return lit_Undef;

    bool sign = false;
    switch (polarity_mode){
    case polarity_true:  sign = false; break;
    case polarity_false: sign = true;  break;
    case polarity_user:  sign = polarity[next]; break;
    case polarity_rnd:   sign = irand(random_seed, 2); break;
    default: assert(false); }

    if( next != var_Undef && phase_saving && phase[next] != l_Undef )
      sign = (phase[next] == l_False);

    return next == var_Undef ? lit_Undef : Lit(next, sign);
}

Lit Solver::pickBranchLitLex()
{
  for(int i = 0; i != cspvars.size(); ++i) {
    cspvar x(i);
    if( x.min(*this) != x.max(*this) ) {
      return pickBranchLitFrom(x);
    }
  }
  // everything instantiated. use pickBranchLitVSIDS to finish off any
  // propositional strugglers
  return pickBranchLitVSIDS(polarity_mode, random_var_freq);
}

Lit Solver::pickBranchLitDom()
{
  int mind{0};
  int minv{-1};
  for(int i = 0; i != cspvars.size(); ++i) {
    cspvar x(i);
    if(x.min(*this) != x.max(*this) && (minv == -1 || x.domsize(*this) < mind)) {
      minv = i;
      mind = x.domsize(*this);
    }
  }
  if (minv != -1)
    return pickBranchLitFrom(cspvar{minv});
  return pickBranchLitVSIDS(polarity_mode, random_var_freq);
}

Lit Solver::pickBranchLitFrom(cspvar x)
{
  switch(valbranch) {
  case VAL_VSIDS: assert(0);
  case VAL_LEX:   return Lit( x.eqi(*this, x.min(*this)) );
  case VAL_BISECT:assert(0);
  }
  assert(0);
}

Lit Solver::pickBranchLit(int polarity_mode, double random_var_freq)
{
  switch(varbranch) {
  case VAR_VSIDS: return pickBranchLitVSIDS(polarity_mode, random_var_freq);
  case VAR_LEX:   return pickBranchLitLex();
  case VAR_DOM:   return pickBranchLitDom();
  case VAR_DOMWDEG:
    assert(0);
  case VAR_USER: {
      assert(user_brancher);
      user_candidates.clear();
      user_brancher(user_candidates);
      if (user_candidates.empty())
        return pickBranchLitVSIDS(polarity_mode, random_var_freq);
      return *std::max_element(
          begin(user_candidates), end(user_candidates),
          [&](Lit a, Lit b) { return activity[var(a)] < activity[var(b)]; });
  } break;
  }
  assert(0);
}

lbool Solver::currentVarPhase(Var x) const { return phase[x]; }

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|
|  Description:
|    Analyze conflict and produce a reason clause.
|
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|
|  Effect:
|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
|________________________________________________________________________________________________@*/
void Solver::analyze(Clause* inconfl, vec<Lit>& out_learnt, int& out_btlevel)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // if a propagator is not monotone (i.e., it may happen f(D1)
    // \subset f(D2) even though D2 \subset D1), it may fail with a
    // clause that contains no literals at the current decision level,
    // leaving poor analyze all confused. This may even mean that the
    // conflict clause is empty.
    int maxlvl = 0;
    for (int i = 0; i != inconfl->size(); ++i) {
        Lit l = (*inconfl)[i];
        assert(var(l) != var_Undef);
        assert(value(l) == l_False);
        if (level[var(l)] > maxlvl)
            maxlvl = level[var(l)];
    }
    if (trace && maxlvl < decisionLevel()) {
      cout << "clause failed at level " << maxlvl
           << ", backtracking before analysis\n";
    }
    cancelUntil(maxlvl);
    if (maxlvl == 0) {
        out_btlevel = 0;
        return;
    }

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;
    out_btlevel = 0;

    if( trace && debugclauses )
      cout << "cnfl analysis starting with "
           << print(*this, inconfl) << "\n";

    explanation_ptr confl{inconfl};

    do{
        assert(!confl.null());          // (otherwise should be UIP)
        auto resolve_with = [&](auto &expl) {
          if (trace && debugclauses)
              cout << "resolving on " << lit_printer(*this, p) << " with "
                   << print(*this, &expl) << "\n";

          for (Lit q : expl) {
            if (p == q)
              continue;

            if (!seen[var(q)] && level[var(q)] > 0) {
              varBumpActivity(var(q));
              seen[var(q)] = 1;
              if (level[var(q)] >= decisionLevel())
                pathC++;
              else {
                out_learnt.push(q);
                if (level[var(q)] > out_btlevel)
                  out_btlevel = level[var(q)];
              }
            }
          }
        };

        if (confl.has<Clause>()) {
          auto &c = *confl.get<Clause>();

          if (c.learnt())
            claBumpActivity(c);

          resolve_with(c);
        } else {
            auto *e = confl.get<explainer>();
            auto &explbuffer = analyze_explbuffer;
            explbuffer.clear();
            if (trace && debugclauses)
                cout << "generating deferred clause from " << *e << " for "
                     << lit_printer(*this, p) << "\n";
            e->explain(*this, p, explbuffer);
            assert(explbuffer[0] == p);

            resolve_with(explbuffer);
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason[var(p)];
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    if( trace && debugclauses )
      cout << "before minimization "
           << print(*this, &out_learnt) << "\n";

    // Simplify conflict clause:
    //
    int i, j;
    if (expensive_ccmin){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        out_learnt.copyTo(analyze_toclear);
        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason[var(out_learnt[i])].null() || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
    }else{
        out_learnt.copyTo(analyze_toclear);
        for (i = j = 1; i < out_learnt.size(); i++){
            Clause& c = *explicit_reason(~out_learnt[i]);
            for (int k = 1; k < c.size(); k++)
                if (!seen[var(c[k])] && level[var(c[k])] > 0){
                    out_learnt[j++] = out_learnt[i];
                    break; }
        }
    }
    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    reduce_var(out_learnt);
    tot_literals += out_learnt.size();

    if( trace && debugclauses )
      cout << "after minimization "
           << print(*this, &out_learnt) << "\n";

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        for (int i = 2; i < out_learnt.size(); i++)
            if (level[var(out_learnt[i])] > level[var(out_learnt[max_i])])
                max_i = i;
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level[var(p)];
    }


    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
#ifdef INVARIANTS
    for(int i = 0; i != nVars(); ++i) assert(seen[i] == 0);
#endif
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason[var(analyze_stack.last())]);
        Clause& c = *explicit_reason(~analyze_stack.last());
        analyze_stack.pop();

        for (int i = 0; i < c.size(); i++){
            if( c[i] == p ) continue;
            Lit p  = c[i];
            if (!seen[var(p)] && level[var(p)] > 0){
                if (!reason[var(p)].null() && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}

// If a literal has an explicit Clause * as its reason, return
// this. Otherwise, make it explicit and return it
Clause *Solver::explicit_reason(Lit p) {
  auto& r = reason[var(p)];
  if (r.null())
      return nullptr;
  if (r.has<Clause>())
    return r.get<Clause>();

  auto *e = r.get<explainer>();
  auto &explbuffer = analyze_explbuffer;
  explbuffer.clear();
  e->explain(*this, p, explbuffer);
  assert(explbuffer[0] == p);

  if (trace && debugclauses)
      cout << "generated deferred clause " << print(*this, &explbuffer)
           << " from " << *e << " for " << lit_printer(*this, p) << "\n";

  if (debugclauses)
      for(Lit q : explbuffer)
          assert(q == p
              || (value(q) == l_False && level[var(q)] <= level[var(p)]));

  auto *cl = Clause_new(explbuffer);
  r.set(cl);
  addInactiveClause(cl);
  return cl;
}

// make sure a clause has no redundant information wrt to any one
// variable: if it contains Xi != d, it cannot contain any other
// literal of Xi. If it contains Xi <= d, it cannot containt Xi=d' or
// Xi <= d' with d' < d. If it contains ~(Xi <= d), it cannot contain
// Xi=d' or ~(Xi <= d') with d' > d.
//
// This assumes the clause is not contradictory, i.e., does not
// contain both Xi != d and Xi != d' for d != d'
//
// This is guaranteed to not remove the UIP. Proof: it only removes
// literals that are implied by other literals that are already in the
// clause. If the UIP is removed, it means that it is implied by a
// literal l. But since l implies the UIP, it must be in the same
// decision level, contradiction.
void Solver::reduce_var(vec<Lit>& ps)
{
  vector<char> & seen = reduce_var_seen;
  vector<int> & varmin = reduce_var_min;
  vector<int> & varmax = reduce_var_max;
  vector<char> & varasgn = reduce_var_asgn;
  vector<cspvar> & toclear = reduce_var_toclear;

  // note here that literals are false, so the literal with event LEQ
  // is false when the var is > d (>= d+1)
  for(int i = 0; i != ps.size(); ++i) {
    domevent de = events[toInt(ps[i])];
    if( de.type != domevent::NONE )
      if( !seen[de.x._id] ) {
        seen[de.x._id] = true;
        toclear.push_back(de.x);
      }
    switch( de.type ) {
    case domevent::NONE:
    case domevent::EQ:
      break;
    case domevent::NEQ:
      varmin[ de.x._id ] = de.d;
      varmax[ de.x._id ] = de.d;
      varasgn[de.x._id ] = true;
      break;
    case domevent::GEQ:
      varmax[de.x._id ] = min(varmax[de.x._id], de.d-1);
      break;
    case domevent::LEQ:
      varmin[de.x._id ] = max(varmin[de.x._id], de.d+1);
      break;
    }
  }

  int i = 0, j = 0;
  for(; i != ps.size(); ++i) {
    domevent de = events[toInt(ps[i])];
    switch( de.type ) {
    case domevent::EQ:
      if( de.d >= varmin[de.x._id] && de.d <= varmax[de.x._id] )
        ps[j++] = ps[i];
      break;
    case domevent::NONE:
      ps[j++] = ps[i];
      break;
    case domevent::NEQ:
      ps[j++] = ps[i];
      break;
    case domevent::GEQ:
      if( !varasgn[de.x._id] && de.d == varmax[de.x._id]+1 )
        ps[j++] = ps[i];
      break;
    case domevent::LEQ:
      if( !varasgn[de.x._id] && de.d == varmin[de.x._id]-1 )
        ps[j++] = ps[i];
      break;
    }
  }
  ps.shrink(i-j);

  for(size_t i = 0; i != toclear.size(); ++i) {
    cspvar x = toclear[i];
    varmin[ x._id ] = cspvars[x._id].omin;
    varmax[ x._id ] = cspvars[x._id].omax;
    varasgn[ x._id ] = false;
    seen[ x._id] = false;
  }
  toclear.clear();
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason[x].null()){
                assert(level[x] > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = *explicit_reason(trail[i]);
                for (int j = 1; j < c.size(); j++)
                    if (level[var(c[j])] > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}

void Solver::debugclause(Clause *from, cons *c)
{
  if( !allow_clause_dbg ) return;

  if(trace) {
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    cout << "Debugging clause " << print(*this, from)
         << " generated by constraint "
         << cons_state_printer(*this, *c) << "\n";
  }

  Solver s1;
  s1.allow_clause_dbg = false; // to avoid infinite recursion
  // add all variables
  int nv = nVars();
  s1.watches.growTo(2*nv);
  s1.wakes_on_lit.growTo(nv);
  s1.sched_on_lit.growTo(nv);
  s1.reason.growTo(nv);
  s1.assigns.growTo(nv, toInt(l_Undef));
  s1.level.growTo(nv, -1);
  s1.activity.growTo(nv, 0);
  s1.seen.growTo(nv, 0);

  polarity.copyTo(s1.polarity);
  decision_var.copyTo(s1.decision_var);

  s1.phase.growTo(nv, l_Undef);
  events.copyTo(s1.events);
  setevents.copyTo(s1.setevents);

  s1.varnames.resize(nv);
  s1.cspvarnames.resize( cspvarnames.size() );
  s1.setvarnames.resize( setvarnames.size() );

  cspvars.copyTo(s1.cspvars);
  for(int i = 0; i != cspvars.size(); ++i) {
    cspvar x(i);
    cspvar_fixed & xf = s1.cspvars[i];
    xf.min = xf.omin;
    xf.max = xf.omax;
    xf.dsize = xf.max - xf.min + 1;
    if( xf.min == xf.max ) {
      s1.uncheckedEnqueue_np(Lit(cspvars[i].firstbool+1), NO_REASON);
      s1.uncheckedEnqueue_np(Lit(cspvars[i].firstbool), NO_REASON);
    } else
      x.setmax(s1, xf.max, (Clause*)0L);
  }

  setvars.copyTo(s1.setvars);

  s1.propagate();

  s1.trace = trace;
  c->clone(s1);
  for(int i = 0; i != from->size(); ++i) {
    Lit q = (*from)[i];
    if(s1.value(q) == l_True )
      return; // conflict already
    if( s1.value( q ) != l_False )
      s1.uncheckedEnqueue(~q, 0L);
  }
  Clause *confl = s1.propagate();
  if(trace) {
    cout << "clause debugging: constraint is now in state\n"
         << cons_state_printer(s1, *s1.conses[0]) << "\n";
    if( !confl )
      cout << "*** no failure ***\n";
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  }
  if( !confl ) {
    cout << "constraint " << cons_state_printer(*this, *c)
         << "\ngenerated incorrect clause\n";
  }
  assert(confl);
}

void Solver::check_debug_solution(Lit p, explanation_ptr from)
{
#ifdef INVARIANTS
    assert(!active_constraint || !learning || from);
    // decisions are allowed to be wrong
    if (learning && from.null() && decisionLevel() > 0)
        return;
    // when not learning, a decision is the first literal at the current dlvl
    if (!learning && decisionLevel() > 0
        && trail.size() == trail_lim[decisionLevel()-1])
        return;

    if (!debug_solution_lits.empty() &&
        (p == lit_Undef ||
         sign(p) != (debug_solution_lits[var(p)] == l_False))) {
        // p disagrees with stored solution
        bool consistent{true};
        for(int i = 0; i != nVars(); ++i) {
            if (i != var(p) && value(i) != l_Undef &&
                value(i) != debug_solution_lits[i]) {
                consistent = false;
                break;
            }
        }
        if (consistent)
            cout << "p = " << lit_printer(*this, p)
                 << " is inconsistent, rest of assignment consistent = "
                 << consistent << " reason = " << from
                 << " dlvl = " << decisionLevel() << endl;
        assert(!consistent);
    }
    if( debug_solution.empty() ) return;
    domevent pe;
    if( p != lit_Undef ) {
        pe = events[toInt(p)];
        switch(pe.type) {
        case domevent::NEQ:
            if( pe.d != debug_solution[pe.x._id] ) return;
            break;
        case domevent::EQ:
            if( pe.d == debug_solution[pe.x._id] ) return;
            break;
        case domevent::LEQ:
            if( pe.d >= debug_solution[pe.x._id] ) return;
            break;
        case domevent::GEQ:
            if( pe.d <= debug_solution[pe.x._id] ) return;
            break;
        case domevent::NONE:
            return;
        }
    }
    cout << "\t\tinconsistent with x"
         << pe.x._id << " == " << debug_solution[pe.x._id] << "\n";
    for(int i = 0; i != cspvars.size(); ++i) {
      if( i != pe.x._id &&
          value(cspvars[i].eqi( debug_solution[i] )) == l_False ) {
        cout << "\t\tx" << i << " != " << debug_solution[i] << " already "
             << "at level " << level[cspvars[i].eqi( debug_solution[i] )]
             << "\n";
        cout << "\t\t\tother assignments inconsistent with solution\n";
        return;
      }
    }
    assert(false);
#endif
}

void Solver::uncheckedEnqueue_np(Lit p, explanation_ptr from)
{
    assert(value(p) == l_Undef);
    assigns [var(p)] = toInt(lbool(!sign(p)));  // <<== abstract but not uttermost effecient
    level   [var(p)] = decisionLevel();
    if (reason[var(p)] && reason[var(p)].has<explainer>())
        reason[var(p)].get<explainer>()->release();
    reason  [var(p)] = from;
    trail.push(p);

#ifdef INVARIANTS
    if (from && from.has<Clause>()) {
      Clause& c = *from.get<Clause>();
      bool foundp = false;
      for(int i = 0; i != c.size(); ++i) {
        assert(c[i] != lit_Undef);
        if( c[i] != p )
          assert( value(c[i]) == l_False );
        else
          foundp = true;
      }
      assert(foundp);
    }
#endif
}

/* We do manual domain updates here, rather than wait for
   propagate(). This ensures that the encoding of the domain is always
   immediately consistent with what propagators ask. For example, if a
   propagator does x.setmin(5), the literals (x <= 4), (x <= 3) etc,
   and (x = 4), (x = 3), etc will all be false immediately with the
   correct reasons.

   A minor advantage of this is that it's probably faster than doing
   unit propagation and does not bring the clauses ps1...ps4 into the
   cache.
 */
void Solver::uncheckedEnqueue_common(Lit p, explanation_ptr from)
{
    if (trace) {
        domevent const& pevent = events[toInt(p)];
        if (!noevent(pevent))
            cout << domevent_printer(*this, pevent);
        else
            cout << lit_printer(*this, p);
        if (active_constraint) {
            cout << " forced by "
                 << cons_state_printer(*this, *active_constraint);
        }
        if (from) {
            cout << " forced by ";
            if (from.has<Clause>())
                cout << " clause " << print(*this, from.get<Clause>());
            else
                cout << " " << *from.get<explainer>();
        }
        cout << " at level " << decisionLevel() << "\n";
    }

    check_debug_solution(p, from);
    uncheckedEnqueue_np(p, from);

#ifdef INVARIANTS
    if (debugclauses && active_constraint && from.has<Clause>())
        debugclause(from.get<Clause>(), active_constraint);
#endif

    // update csp var and propagate, if applicable
    domevent const &pevent = events[toInt(p)];
    if( noevent(pevent) ) return;
    cspvar_fixed& xf = cspvars[pevent.x._id];
    if( pevent.type == domevent::EQ ) {
      xf.dsize = 1;
      xf.min = xf.max = pevent.d;
      // propagate towards max
      if( value(xf.leqiUnsafe(pevent.d)) != l_True ) {
        uncheckedEnqueue_np( Lit(xf.leqiUnsafe(pevent.d)),
                             xf.ps2[pevent.d-xf.omin] );
        if( value(xf.eqiUnsafe(pevent.d+1)) != l_False)
          uncheckedEnqueue_np( ~Lit(xf.eqiUnsafe(pevent.d+1)),
                               xf.ps3[pevent.d+1-xf.omin] );
        int leq = pevent.d+1;
        while( leq < xf.omax && value(xf.leqiUnsafe(leq)) != l_True ) {
          uncheckedEnqueue_np( Lit(xf.leqiUnsafe(leq)), xf.ps1[leq-1-xf.omin] );
          if( value(xf.eqiUnsafe(leq+1)) != l_False)
            uncheckedEnqueue_np( ~Lit(xf.eqiUnsafe(leq+1)), xf.ps3[leq+1-xf.omin] );
          ++leq;
        }
      }
      // propagate towards min
      if( pevent.d > xf.omin &&
          value(xf.leqiUnsafe(pevent.d-1)) != l_False ) {
        uncheckedEnqueue_np( ~Lit(xf.leqiUnsafe(pevent.d-1)),
                             xf.ps3[pevent.d-xf.omin] );
        if( value(xf.eqiUnsafe(pevent.d-1)) != l_False)
          uncheckedEnqueue_np( ~Lit(xf.eqiUnsafe(pevent.d-1)),
                               xf.ps2[pevent.d-1-xf.omin] );
        int geq = pevent.d-2;
        while( geq >= xf.omin && value(xf.leqiUnsafe(geq)) != l_False ) {
          uncheckedEnqueue_np( ~Lit(xf.leqiUnsafe(geq)), xf.ps1[geq-xf.omin] );
          if( value(xf.eqiUnsafe(geq)) != l_False )
            uncheckedEnqueue_np( ~Lit(xf.eqiUnsafe(geq)), xf.ps2[geq-xf.omin] );
          --geq;
        }
      }
    }

    if( pevent.type == domevent::NEQ ) {
      --xf.dsize;
      if( pevent.d == xf.min ) {
        while( value(xf.eqiUnsafe(xf.min)) == l_False ) {
          assert( value(xf.leqiUnsafe(xf.min)) != l_False &&
                  (xf.min == xf.omin ||
                   value(xf.leqiUnsafe(xf.min-1)) == l_False ));
          uncheckedEnqueue_np( ~Lit( xf.leqiUnsafe(xf.min) ),
                               xf.ps4[xf.min-xf.omin] );
          ++xf.min;
        }
      }
      if( pevent.d == xf.max ) {
        while( value(xf.eqiUnsafe(xf.max)) == l_False ) {
          assert( (xf.max == xf.omax ||
                   value(xf.leqiUnsafe(xf.max)) == l_True) &&
                  value(xf.leqiUnsafe(xf.max-1)) != l_False );
          uncheckedEnqueue_np( Lit( xf.leqiUnsafe(xf.max-1) ),
                               xf.ps4[xf.max-xf.omin] );
          --xf.max;
        }
      }
    }

    if( pevent.type == domevent::GEQ ) {
      int geq = pevent.d-1;
      if( value(xf.eqiUnsafe(geq)) != l_False ) {
        uncheckedEnqueue_np( ~Lit( xf.eqiUnsafe(geq)), xf.ps2[geq-xf.omin] );
        --xf.dsize;
      }
      while(geq > xf.omin &&
            value( xf.leqiUnsafe(geq-1) ) != l_False ) {
        uncheckedEnqueue_np( ~Lit(xf.leqiUnsafe(geq-1)),
                             xf.ps1[geq-1-xf.omin] );
        if( value(xf.eqiUnsafe(geq-1)) != l_False ) {
          uncheckedEnqueue_np( ~Lit( xf.eqiUnsafe(geq-1)), xf.ps2[geq-1-xf.omin] );
          --xf.dsize;
        }
        --geq;
      }
      xf.min = pevent.d;
      while( value(xf.eqiUnsafe(xf.min)) == l_False ) {
        assert( value(xf.leqiUnsafe(xf.min)) != l_False &&
                value(xf.leqiUnsafe(xf.min-1)) == l_False );
        uncheckedEnqueue_np( ~Lit( xf.leqiUnsafe(xf.min) ),
                             xf.ps4[xf.min-xf.omin] );
        ++xf.min;
      }
    }
    if( pevent.type == domevent::LEQ ) {
      int leq = pevent.d+1;
      if( leq <= xf.omax && value(xf.eqiUnsafe(leq)) != l_False ) {
        uncheckedEnqueue_np( ~Lit( xf.eqiUnsafe(leq)), xf.ps3[leq-xf.omin] );
        --xf.dsize;
      }
      while( leq < xf.omax &&
             value( xf.leqiUnsafe(leq) ) != l_True ) {
        uncheckedEnqueue_np( Lit( xf.leqiUnsafe(leq) ), xf.ps1[leq-1-xf.omin] );
        if( value(xf.eqiUnsafe(leq+1)) != l_False ) {
          uncheckedEnqueue_np( ~Lit( xf.eqiUnsafe(leq+1)), xf.ps3[leq+1-xf.omin] );
          --xf.dsize;
        }
        ++leq;
      }
      xf.max = pevent.d;
      while( value(xf.eqiUnsafe(xf.max)) == l_False ) {
        assert( value(xf.leqiUnsafe(xf.max)) == l_True &&
                value(xf.leqiUnsafe(xf.max-1)) != l_False );
        uncheckedEnqueue_np( Lit( xf.leqiUnsafe(xf.max-1) ),
                             xf.ps4[xf.max-xf.omin] );
        --xf.max;
      }
    }

    if( xf.max == xf.min &&
        value( xf.eqiUnsafe(xf.max) ) != l_True )
      uncheckedEnqueue_np( Lit(xf.eqiUnsafe(xf.max)),
                           xf.ps4[xf.max-xf.omin] );

#ifdef INVARIANTS
    for(int i = xf.omin; i != xf.omax; ++i ) {
      if( value(xf.leqiUnsafe(i)) == l_False )
        assert( xf.min > i );
      if( value(xf.leqiUnsafe(i)) == l_True )
        assert( xf.max <= i );
      if( value(xf.eqiUnsafe(i)) == l_False )
        assert( xf.min != i && xf.max != i );
    }
#endif
}

void Solver::uncheckedEnqueue(Lit p, Clause* from)
{
    uncheckedEnqueue_common(p, explanation_ptr(from));
}

void Solver::uncheckedEnqueueDeferred(Lit p, explainer *from)
{
    uncheckedEnqueue_common(p, explanation_ptr(from));
    from->use();
}

Clause *Solver::enqueueFill(Lit p, vec<Lit>& ps)
{
  if( value(p) == l_True ) return 0L;
  if (debugclauses)
    for (int i = 0; i != ps.size(); ++i) {
      if (value(ps[i]) != l_False) {
        cout << "Trying to force " << lit_printer(*this, p)
             << " with invalid clause " << print(*this, &ps) << endl;
        assert(value(ps[i]) == l_False);
      }
    }
  ps.push(p);
  Clause *r = Clause_new(ps, false, p);
  ps.pop();
  addInactiveClause(r);
  if( value(p) == l_False ) return r;
  uncheckedEnqueue( p, r );
  return 0L;
}

void Solver::nonMonotoneEnqueue(Lit p, Clause *from)
{
    assert(value(p) == l_Undef);
    assert(from);
    assert((*from)[0] == p);
    if (trace) {
      int maxlevel =
          from->size() == 1
              ? 0
              : varLevel(*std::max_element(
                    begin(*from) + 1, end(*from), [&](auto l1, auto l2) {
                      return this->varLevel(l1) < this->varLevel(l2);
                    }));
      cout << "Backpruning " << lit_printer(*this, p) << " to level "
           << maxlevel << " with clause " << print(*this, from) << "\n";
    }
    newDecisionLevel();
    uncheckedEnqueue(~p, NO_REASON);
}

/*************************************************************************

  Queue handling

  schedule(c): puts the constraint c into the queue it wants. Note
  this is static, otherwise we would need a virtual function call +
  bringing the cons into the cache

  unschedule(c): remove c from the queue

  reset_queue(): remove all constraints from the queue
 */
void Solver::schedule(int cqidx)
{
  consqueue & cq = consqs[cqidx];
  consqueue & p = consqs[cq.priority+1];
  if( cq.next >= 0 ) return; // already scheduled
  cq.prev = p.prev;
  cq.next = cq.priority+1;
  p.prev = cqidx;
  consqs[cq.prev].next = cqidx;
}

int Solver::first_scheduled()
{
  int cqidx = 0;
  while( cqidx >= 0 && !consqs[cqidx].c )
    cqidx = consqs[cqidx].next;
  if( cqidx >= 0 && consqs[cqidx].c ) return cqidx;
  else return -1;
}

void Solver::unschedule(int cqidx)
{
  consqueue * cq = &consqs[cqidx];
  consqs[cq->next].prev = cq->prev;
  consqs[cq->prev].next = cq->next;
  cq->next = -1;
  cq->prev = -1;
}

void Solver::reset_queue()
{
  int cqidx = 0;
  while( cqidx >= 0 ) {
    int nxt = consqs[cqidx].next;
    consqs[cqidx].next = -1;
    consqs[cqidx].prev = -1;
    cqidx = nxt;
  }
  for(int i = 0; i <= MAX_PRIORITY+1; ++i) {
    if( i <= MAX_PRIORITY )
      consqs[i].next = i+1;
    if( i )
      consqs[i].prev = i-1;
  }
}

Clause *Solver::propagate_wakes(Lit p, const vec<wake_stub>& wakes)
{
  Clause *confl = NO_REASON;
  /* Now propagate constraints that wake on this literal */
  for (const wake_stub &ws : wakes) {
    cons *con = ws.first;
    active_constraint = con;
    if (ws.second)
      confl = con->wake_advised(*this, p, ws.second);
    else
      confl = con->wake(*this, p);
    active_constraint = 0L;
    if (confl) {
      if (trace) {
        cout << "Constraint " << cons_state_printer(*this, *con) << " failed, "
             << "clause " << print(*this, confl) << endl;
      }
      if (debugclauses)
        debugclause(confl, con);
      qhead = trail.size();
      break;
    }
  }
  return confl;
}

/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise NULL.
|
|    Propagation involves both unit propagation and constraint propagation
|    for constraints that wake on changes.
|
|    Additionally, it puts in the constraint queue all propagators that
|    want to be scheduled on changes.
|
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
Clause* Solver::propagate_inner()
{
    Clause* confl     = NULL;
    int     num_props = 0;

    while (qhead < trail.size()){
        Lit p   = trail[qhead++];     // 'p' is enqueued fact to propagate.

        // binary clauses first
        vec<watch>& bws = binwatches[toInt(p)];
        for(watch w : bws) {
            Lit other = w.block;
            assert(other != lit_Undef);
            if (value(other) == l_False)
                return w.c;
            else if (value(other) != l_True)
                uncheckedEnqueue(other, w.c);
        }


        vec<watch>& ws  = watches[toInt(p)];
        vec<watch>::iterator i, j, end;
        num_props++;

        for (i = j = ws.begin(), end = ws.end(); i != end;) {
          watch &w = *i;

          Lit blocker = w.block;
          if (blocker != lit_Undef && value(blocker) == l_True) {
            *j++ = *i++;
            continue;
          }

          ++i;
          Clause &c = *w.c;

          // Make sure the false literal is data[1]:
          Lit false_lit = ~p;
          if (c[0] == false_lit)
            std::swap(c[0], c[1]);

          assert(c[1] == false_lit);

          // If 0th watch is true, then clause is already satisfied.
          Lit first = c[0];
          if (first != blocker && value(first) == l_True) {
            *j++ = w;
            continue;
          }

          // Look for new watch:
          if (UP_callback) {
              UP_result upr = UP_callback(p, c);
              if (upr.newwatchidx >= 0) {
                  if (upr.newwatchidx <= 1)
                    throw std::logic_error(
                        "New watch cannot be one of the existing watches");
                  std::swap(c[1], c[upr.newwatchidx]);
                  watches[toInt(~c[1])].push(w);
                  goto FoundWatch;
              }
          } else {
            for (int k = 2; k < c.size(); k++)
              if (value(c[k]) != l_False) {
                c[1] = c[k];
                c[k] = false_lit;
                watches[toInt(~c[1])].push(w);
                goto FoundWatch;
              }
          }

          // Did not find watch -- clause is unit under assignment:
          *j++ = w;
          if (value(first) == l_False) {
            if (trace)
              cout << "Clause " << print(*this, &c) << " failed\n";
            confl = w.c;
            qhead = trail.size();
            // Copy the remaining watches:
            while (i < end)
              *j++ = *i++;
          } else
            uncheckedEnqueue(first, &c);
        FoundWatch:;
        }
        ws.shrink(i - j);

        if(confl)
          break;

        confl = propagate_wakes(p, wakes_on_lit[var(p)]);
        if(confl)
          break;

        for (int consid : sched_on_lit[var(p)])
            schedule(consid);

        domevent const & pe = events[toInt(p)];
        vec< pair<cons*, void*> > *dewakes = 0L;
        vec<int> *desched = 0L;
        switch(pe.type) {
        case domevent::NEQ:
          dewakes=&(cspvars[pe.x._id].wake_on_dom);
          desched=&(cspvars[pe.x._id].schedule_on_dom);
          break;
        case domevent::EQ:
          dewakes=&(cspvars[pe.x._id].wake_on_fix);
          desched=&(cspvars[pe.x._id].schedule_on_fix);
          break;
        case domevent::LEQ:
          if( pe.x.max(*this) >= pe.d ) {
            dewakes=&(cspvars[pe.x._id].wake_on_ub);
            desched=&(cspvars[pe.x._id].schedule_on_ub);
          }
          break;
        case domevent::GEQ:
          if( pe.x.min(*this) <= pe.d ) {
            dewakes=&(cspvars[pe.x._id].wake_on_lb);
            desched=&(cspvars[pe.x._id].schedule_on_lb);
          }
          break;
        case domevent::NONE:
          break;
        }

        if (dewakes) {
          confl = propagate_wakes(p, *dewakes);
          if (confl)
            break;

          for (int consid : *desched)
            schedule(consid);
        }

        setevent const &se = setevents[toInt(p)];
        if (noevent(se))
          continue;
        vec<wake_stub> *sewakes = 0L;
        vec<int> *sesched = 0L;
        switch (se.type) {
        case setevent::IN:
          sewakes = &(setvars[se.x._id].wake_on_in);
          sesched = &(setvars[se.x._id].schedule_on_in);
          break;
        case setevent::EX:
          sewakes = &(setvars[se.x._id].wake_on_ex);
          sesched = &(setvars[se.x._id].schedule_on_ex);
          break;
        case setevent::NONE:
          assert(0);
        }

        if (sewakes) {
          confl = propagate_wakes(p, *sewakes);
          if (confl)
            break;

          for (int sconsid : *sesched)
            schedule(sconsid);
        }
    }
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}

Clause *Solver::propagate()
{
  Clause *confl = NULL;
  int next = -1;
  do {
    if( next >= 0 ) {
      unschedule(next);
      active_constraint = consqs[next].c;
      confl = consqs[next].c->propagate(*this);
      active_constraint = 0L;
      if( confl ) {
        if( trace ) {
          cout << "Constraint "
               << cons_state_printer(*this, *consqs[next].c) << " failed, "
               << "clause @ " << print(*this, confl) << "\n";
        }
        if( debugclauses )
          debugclause(confl, consqs[next].c);
        qhead = trail.size();
        reset_queue();
        return confl;
      }
    }
    if( qhead < trail.size() )
      confl = propagate_inner();
    if( confl ) {
      reset_queue();
      return confl;
    }

    next = first_scheduled();
  } while( qhead < trail.size() || next >= 0);
  return 0L;
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { bool operator () (Clause* x, Clause* y) { return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity()); } };
void Solver::reduceDB()
{
    int     i, j;
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

    sort(learnts, reduceDB_lt());
    for (i = j = 0; i < learnts.size() / 2; i++){
        if (learnts[i]->size() > 2 && !locked(*learnts[i]))
            removeClause(*learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    for (; i < learnts.size(); i++){
        if (learnts[i]->size() > 2 && !locked(*learnts[i]) && learnts[i]->activity() < extra_lim)
            removeClause(*learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);
}

void Solver::gcInactive()
{
  int i,j;
  for(i = j = 0; i < inactive.size(); ++i) {
    if( !locked(*inactive[i]) )
      free(inactive[i]);
    else
      inactive[j++] = inactive[i];
  }
  inactive.shrink(i - j);
}

void Solver::removeSatisfied(vec<Clause*>& cs)
{
    int i,j;
    for (i = j = 0; i < cs.size(); i++){
        if (satisfied(*cs[i]))
            removeClause(*cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != NULL)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);

    // Remove fixed variables from the variable heap:
    order_heap.filter(VarFilter(*this));

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]
|
|  Description:
|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
|    indicate infinity.
|
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts, double* nof_learnts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;

    starts++;

    for (;;){
        Clause *confl;
        try {
          confl = propagate();
        } catch( unsat& ) {
          assert(decisionLevel() == 0);
          ++conflicts;
          if (trace)
              cout << "UNSAT after exception\n";
          return l_False;
        }
        if (confl != NULL){
            // CONFLICT
            check_debug_solution(lit_Undef, confl);
            conflicts++; conflictC++;
            if (decisionLevel() == 0) {
                if (trace)
                    cout << "UNSAT after conflict at dlvl 0\n";
                return l_False;
            }

            if (!learning) {
              if (!clause_callbacks.empty()) {
                  // we learn no clause, only implicitly the set of
                  // all decisions
                  learnt_clause.clear();
                  learnt_clause.push(~decisionAtLevel(decisionLevel()));
                  for (int i = 1; i < decisionLevel(); ++i)
                      learnt_clause.push(~decisionAtLevel(i));
                  for (auto& f : clause_callbacks) {
                      auto r = f(learnt_clause, decisionLevel() - 1);
                      if (r != CCB_OK)
                        interrupt_requested = true;
                      if (r == CCB_MODIFIED)
                        cout << "Warning: clause modified but learning not "
                                "enabled\n";
                  }
                  if (interrupt_requested || !withinBudget()) {
                      cancelUntil(0);
                      return l_Undef;
                  }
              }
              Lit flip = trail[ trail_lim[ decisionLevel() - 1 ] ];
              cancelUntil( decisionLevel() - 1 );
              uncheckedEnqueue(~flip, 0L);
              continue;
            }

            bool reanalyze{false};
            bool non_asserting{false};
            do {
              reanalyze = false;
              learnt_clause.clear();
              analyze(confl, learnt_clause, backtrack_level);

              for (auto &f : clause_callbacks) {
                auto r = f(learnt_clause, backtrack_level);
                switch (r) {
                case CCB_OK:
                  break;
                case CCB_INTERRUPT:
                  interrupt_requested = true;
                  break;
                case CCB_MODIFIED:
                  for (Lit l : learnt_clause)
                    assert(value(l) == l_False);
                  assert(
                      std::none_of(begin(learnt_clause), end(learnt_clause),
                                   [&](Lit l) { return value(l) == l_True; }));
                  erase_if(learnt_clause, [&](Lit l) {
                    return (value(l) != l_Undef && varLevel(l) == 0);
                  });

                  if (learnt_clause.size() == 0)
                    return l_False;
                  if (learnt_clause.size() == 1)
                    backtrack_level = 0;
                  else {
                    std::partial_sort(begin(learnt_clause),
                                      begin(learnt_clause) + 2,
                                      end(learnt_clause), [&](Lit a, Lit b) {
                                        return varLevel(a) > varLevel(b);
                                      });
                    if (backtrack_level != varLevel(learnt_clause[1]))
                      backtrack_level =
                          std::min(backtrack_level, varLevel(learnt_clause[1]));
                    cancelUntil(backtrack_level);
                    if (value(learnt_clause[0]) == l_False &&
                        varLevel(learnt_clause[0]) ==
                            varLevel(learnt_clause[1])) {
                      confl = addInactiveClause(learnt_clause);
                      reanalyze = true;
                    } else if (value(learnt_clause[1]) != l_False)
                        non_asserting = true;
                  }
                  break;
                }
              }
            } while (reanalyze);

            if (interrupt_requested || !withinBudget()) {
              cancelUntil(0);
              return l_Undef;
            }

            cancelUntil(backtrack_level);

            if (learnt_clause.size() == 0) {
              if (trace)
                cout << "UNSAT with empty explanation clause\n";
              return l_False;
            }
            assert(value(learnt_clause[0]) == l_Undef);

            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                Clause* c = Clause_new(learnt_clause, true);
                learnts.push(c);
                attachClause(*c);
                claBumpActivity(*c);

                // a user provided rewriting function can produce
                // non-asserting clauses
                assert(non_asserting || value(learnt_clause[1]) == l_False);
                if (value(learnt_clause[1]) == l_False)
                    uncheckedEnqueue(learnt_clause[0], c);
            }

            if(trace) {
              cout << "Conflict " << conflicts
                   << ", backjumping to level " << backtrack_level
                   << " and setting " << lit_printer(*this, learnt_clause[0]);
              if( !noevent(event(learnt_clause[0])) ) {
                cspvar x = event(learnt_clause[0]).x;
                cout << ", " << cspvar_printer(*this, x) << " in "
                     << domain_as_set(*this, x);
              }
              cout << "\n";
              auto latest_clause = reason[var(learnt_clause[0])];
              assert(latest_clause.has<Clause>());
              auto *cl = latest_clause.get<Clause>();
              cout << "Conflict clause " << conflicts
                   << ": " << print(*this, cl)
                   << "\n";
            }

            varDecayActivity();
            claDecayActivity();

        }else{
            // NO CONFLICT

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if ( decisionLevel() == 0 || inactive.size() > 10*nAssigns() )
              gcInactive();

            /* Increase decision level. We must now call
               newDecisionLevel() immediately after propagate() in
               order to store the changes to backtrackable data made
               by propagation. Otherwise, we might restart next and
               cancelUntil(0) will overwrite our changes. This only
               happens if we restart immediately after learning a unit
               clause.

               This is a new requirement compared to stock minisat
             */
            newDecisionLevel();


            if (restarting &&
                nof_conflicts >= 0 && conflictC >= nof_conflicts){
                // Reached bound on number of conflicts:
                if (trace)
                    cout << "Restarting\n";
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }


            if (*nof_learnts >= 0 && learnts.size()-nAssigns() >= *nof_learnts) {
                // Reduce the set of learnt clauses:
                reduceDB();
                *nof_learnts   *= learntsize_inc;
            }

            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                next = pickBranchLit(polarity_mode, random_var_freq);

                if (next == lit_Undef) {
                    // Model found:
                    if (trace)
                        cout << "Solution found\n";
                    return l_True;
                }
            }

            if( trace ) {
              domevent pe = events[toInt(next)];
              cout << "Decision " << decisions << ":";
              if( noevent(pe) ) {
                cout << lit_printer(*this, next);
              } else {
                cout << domevent_printer(*this, pe);
              }
              cout << " at level " << decisionLevel() << "\n";
            }

            // enqueue 'next'
            assert(value(next) == l_Undef);
            uncheckedEnqueue(next);
            for (auto& f : decision_callbacks)
                f();
        }
    }
}


double Solver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

bool Solver::solve(const vec<Lit>& assumps)
{
    conflict_lim = -1;
    // uninterruptible, because we do not know how to report interruption.
    lbool status{l_Undef};
    while ((status = solveBudget(assumps)) == l_Undef) {
        cerr << "c WARNING: a callback requested interrupt but solver was not "
                "invoked with solveBudget(). Continuing...\n";
        interrupt_requested = false;
    }
    return (status == l_True);
}

lbool Solver::solveBudget(const vec<Lit>& assumps)
{
    model.clear();
    conflict.clear();

    if (!ok) return false;

    interrupt_requested = false;
    assumps.copyTo(assumptions);

    double  nof_conflicts = restart_first;
    double  nof_learnts   = std::max(100.0, nClauses() * learntsize_factor);
    lbool   status        = l_Undef;

    if (verbosity >= 1){
        reportf("============================[ Search Statistics ]==============================\n");
        reportf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        reportf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        reportf("===============================================================================\n");
    }

    int lubycounter = 0;
    int lubybits = 0;
    int lubymult = 1;

    // Search:
    while (status == l_Undef && !interrupt_requested && withinBudget()) {
        if (verbosity >= 1)
            reportf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", (int)conflicts, order_heap.size(), nClauses(), (int)clauses_literals, (int)nof_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progress_estimate*100), fflush(stdout);
        if( !lubybits ) {
          ++lubycounter;
          lubybits = lubycounter ^ (lubycounter-1);
          lubymult = 1;
        }
        //printf("lubymult %d, nof_learnts %f\n", lubymult, nof_learnts);
        nof_conflicts = lubymult * restart_first;
        //nof_conflicts *= restart_inc;
        lubybits >>= 1;
        lubymult <<= 1;
        status = search((int)nof_conflicts, &nof_learnts);
    }

    if (verbosity >= 1)
        reportf("===============================================================================\n");


    if (status == l_True){
        if( trace )
          cout << "Solution ";
        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
        cspmodel.growTo(cspvars.size());
        for(int i = 0; i != cspvars.size(); ++i) {
          cspvar_fixed& xf = cspvars[i];
          cspmodel[i] = std::make_pair(xf.min, xf.max);
          if( trace ) {
            if( i ) cout << ", ";
            cout << cspvar_printer(*this, cspvar(i))
                 << " in " << domain_as_set(*this, cspvar(i));
          }
        }
        if( trace ) cout << "\n";
        cspsetmodel.growTo(setvars.size());
        for(int i = 0; i != setvars.size(); ++i) {
          set<int>& lb = cspsetmodel[i].first;
          set<int>& ub = cspsetmodel[i].second;
          lb.clear(); ub.clear();
          for(int j = setvars[i].min; j <= setvars[i].max; ++j) {
            lbool l = value( setvars[i].ini(j) );
            if( l != l_False )
              ub.insert(j);
            if( l == l_True )
              lb.insert(j);
          }
        }

        if( solution_phase_saving ) {
            for(int i = 0; i != nVars(); ++i)
                phase[i] = toLbool(assigns[i]);
        }
#ifndef NDEBUG
        verifyModel();
#endif
    } else if (status == l_False) {
        if (conflict.size() == 0)
            ok = false;
    }

    cancelUntil(0);
    return status;
}

void Solver::excludeLast()
{
  vec<Lit> exclude;
  for(int i = 0; i != cspvars.size(); ++i) {
    cspvar_fixed & xf = cspvars[i];
    pair<int, int> pi = cspmodel[i];
    assert( pi.first == pi.second );
    exclude.push( ~Lit(xf.eqi( pi.first )) );
  }
  for(int i = 0; i != setvars.size(); ++i) {
    setvar_data & xd = setvars[i];
    pair< set<int>, set<int> > const& pi = cspsetmodel[i];
    set<int> const& lb = pi.first;
    set<int> const& ub = pi.second;
    for(int j = xd.min; j <= xd.max; ++j) {
      if( lb.find(j) != lb.end() )
        exclude.push( ~Lit(xd.ini( j )));
      else if( ub.find(j) == ub.end() )
        exclude.push( Lit(xd.ini(j)));
    }
  }
  if( trace ) {
    cout << "adding exclude clause (";
    for(int i = 0; i != exclude.size(); ++i)
      cout << lit_printer(*this, exclude[i]) << ' ';
    cout << ")\n";
  }
  addClause(exclude);
}

std::vector<Lit> Solver::getImpliedLiterals()
{
    assert(decisionLevel() == 0);
    std::vector<Lit> rv;
    for (int i = 0, e = trail.size(); i != e; ++i)
        rv.push_back(trail[i]);
    return rv;
}

std::optional<std::vector<Lit>> Solver::getImplications(Lit l)
{
    assert(decisionLevel() == 0);
    std::optional<std::vector<Lit>> rv;
    if (value(l) == l_False)
        return rv;
    if (value(l) == l_True) {
        rv = vector<Lit>{};
        return rv;
    }

    newDecisionLevel();
    uncheckedEnqueue(l);
    auto confl = propagate();
    if (confl) {
        cancelUntil(0);
        return rv;
    }

    std::vector<Lit> impl;
    for (int i = trail_lim[0], e = trail.size(); i != e; ++i)
        impl.push_back(trail[i]);
    cancelUntil(0);
    rv = impl;
    return rv;
}

//==================================================
// output

ostream& operator<<(ostream& os, domevent_printer dp)
{
  if( noevent(dp._p) ) {
    os << "noevent";
    return os;
  }

  os << cspvar_printer(dp._s, dp._p.x)
     << ' ' << opstring(dp._p.type) << ' '
     << dp._p.d;
  return os;
}

//=================================================================================================
// Debug methods:


void Solver::verifyModel()
{
    bool failed = false;
    for (int i = 0; i < clauses.size(); i++){
        assert(clauses[i]->mark() == 0);
        Clause& c = *clauses[i];
        for (int j = 0; j < c.size(); j++)
            if (modelValue(c[j]) == l_True)
                goto next;

        reportf("unsatisfied clause: ");
        printClause(*clauses[i]);
        reportf("\n");
        failed = true;
    next:;
    }

    assert(!failed);
}


void Solver::checkLiteralCount()
{
    // Check that sizes are calculated correctly:
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (clauses[i]->mark() == 0)
            cnt += clauses[i]->size();

    if ((int)clauses_literals != cnt){
        fprintf(stderr, "literal count: %d, real value = %d\n", (int)clauses_literals, cnt);
        assert((int)clauses_literals == cnt);
    }
}

Clause* cons::wake_advised(Solver&, Lit, void *)
{
  assert(0);
  return 0L;
}

Clause* cons::wake(Solver&, Lit)
{
  assert(0);
  return 0L;
}

Clause* cons::propagate(Solver&)
{
  assert(0);
  return 0L;
}

void cons::clone(Solver &other)
{
  assert(0);
}

ostream& cons::print(Solver&, ostream& os) const
{
  os << "cons@" << this;
  return os;
}

ostream& cons::printstate(Solver&, ostream& os) const
{
  os << "cons@" << this;
  return os;
}

} //namespace minicsp
