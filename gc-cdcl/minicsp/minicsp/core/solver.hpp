/**-*- Mode:c++; c-basic-offset: 4 -*-******************************************
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

#ifndef __MINICSP_SOLVER_HPP
#define __MINICSP_SOLVER_HPP

#include <vector>
#include <set>
#include <cstdio>
#include <cstring>
#include <string>
#include <functional>
#include <optional>

#include "minicsp/mtl/Vec.h"
#include "minicsp/mtl/Heap.h"
#include "minicsp/mtl/Alg.h"

#include "solvertypes.hpp"

namespace minicsp {

#ifdef _MSC_VER
#define snprintf _snprintf_s
#endif

//=================================================================================================
// Solver -- the main class:

class Solver;
/* a completely opaque reference type. can only be used with Solver::deref */
class btptr {
  size_t offset;
  friend class Solver;
};

/* Branching heuristics: if VSIDS, choose a literal
 * directly. Otherwise, choose a variable, then choose a value using
 * the corresponding value branching heuristic.
 */
enum BranchHeuristic {
    VAR_VSIDS, VAR_LEX, VAR_DOM, VAR_DOMWDEG, VAR_USER
};

enum ValBranchHeuristic {
    VAL_VSIDS, VAL_LEX, VAL_BISECT
};

class Solver {
public:
    struct VarOrderLt {
        const vec<double>&  activity;
        bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
        VarOrderLt(const vec<double>&  act) : activity(act) { }
    };

public:

    // Constructor/Destructor:
    //
    Solver();
    ~Solver();

    // Problem specification:
    //
    Var     newVar    (bool polarity = true, bool dvar = true); // Add a new variable with parameters specifying variable mode.
    cspvar  newCSPVar (int min, int max);                       // Add a CSP (multi-valued) var with the given lower and upper bounds
    std::vector<cspvar> newCSPVarArray(int n, int min, int max);// Add a number of identical CSP vars
    setvar  newSetVar (int min, int max);                       // Add a CSP (multi-valued) var with the given lower and upper bounds
    std::vector<setvar> newSetVarArray(int n, int min, int max);// Add a number of identical CSP vars
    void    addClause (vec<Lit>& ps);                           // Add a clause to the solver. NOTE! 'ps' may be shrunk by this method!
    template<typename veclit>
    void    addClause (veclit&& ps);                           // Add a clause to the solver. NOTE! 'ps' may be shrunk by this method!
    bool    addConstraint(cons *c);                             // Add a constraint. For transfer of ownership only, everything else in the constructor

    /* Waking means that the constraint is called immediately when we
       process the corresponding literal in the assignment
       stack. Scheduling means that when we process the literal, the
       constraint is placed in a queue and then called when the
       assignment stack has been processed. Similar to the difference
       between a demon and a propagator in ilog. Unlike gecode
       advisors, a constraint can prune when it wakes. A constraint
       can also arrange to receive when it wakes a piece of advice.
       If advice == NULL, then the regular, non-advised wake() is
       called, other wake_advised()
    */
    void    wake_on_lit(Var, cons *c, void *advice = 0L);       // Wake this constraint when the Boolear Var is fixed
    void    wake_on_dom(cspvar, cons*c, void *advice = 0L);     // Wake this constraint when a value of cspvar is pruned
    void    wake_on_lb(cspvar, cons*c, void *advice = 0L);      // Wake this constraint when the lb of cspvar is changed
    void    wake_on_ub(cspvar, cons*c, void *advice = 0L);      // Wake this constraint when the ub of cspvar is changed
    void    wake_on_fix(cspvar, cons*c, void *advice = 0L);     // Wake this constraint when the cspvar is assigned

    void    schedule_on_lit(Var, cons*c);                       // Schedule this constraint when the Boolean Var is fixed
    void    schedule_on_dom(cspvar, cons*c);                    // Schedule this constraint when a value of cspvar is pruned
    void    schedule_on_lb(cspvar, cons*c);                     // Schedule this constraint when the lb of cspvar is changed
    void    schedule_on_ub(cspvar, cons*c);                     // Schedule this constraint when the ub of cspvar is changed
    void    schedule_on_fix(cspvar, cons*c);                    // Schedule this constraint when the cspvar is assigned

    btptr   alloc_backtrackable(unsigned size);                 // allocate memory to be automatically restored to its previous contents on backtracking
    void*   get(btptr p);                                       // get direct pointer to  backtrackable mem. for temporary use only
    template<typename T>
    T&      deref(btptr p);                                     // dereference backtrackable mem. for temporary use only
    template<typename T>
    T*      deref_array(btptr p);                               // dereference an array in backtrackable mem. for temporary use only

    // event information
    domevent event(Lit p) const;                                // get the event associated with a literal, if any
    setevent sevent(Lit p) const;                               // get the set variable event associated with a literal

    // variable naming
    // these should all be self-explanatory
    void setVarName(Var v, std::string const& name);
    void setCSPVarName(cspvar v, std::string const& name);
    void setSetVarName(setvar v, std::string const& name);

    std::string const& getVarName(Var v);
    std::string const& getCSPVarName(cspvar v);
    std::string const& getSetVarName(setvar v);

    // Solving:
    //
    bool    simplify     ();                        // Removes already satisfied clauses.
    bool    solve        (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions. Cannot be interrupted
    bool    solve        ();                        // Search without assumptions.
    lbool   solveBudget  (const vec<Lit>& assumps); // Try to solve with assumptions within a budget, return l_Undef if unable
    lbool   solveBudget  ();                        // As above, without assumptions
    bool    okay         () const;                  // FALSE means solver is in a conflicting state
    void    excludeLast  ();                        // add a clause that excludes the last solution

    // Get information from the solver:
    // the set of implied literals at the root
    std::vector<Lit> getImpliedLiterals();
    // the set of all literals newly implied by propagation after
    // setting l to true, meaning this does not include the literals
    // that are already true anyway. Returns an empty optional if
    // setting to true causes a conflict (contrast to an optional
    // containing an empty vector if it has no implications). Note
    // that even if setting l to true causes a conflict, this will not
    // automatically set l to false at the root.
    std::optional<std::vector<Lit>> getImplications(Lit l);

    // Explaining
    // add a clause that will not be propagated, just as a reason or a
    // conflict clause
    void    addInactiveClause(Clause *c);
    // add a clause from a container. Accepts rvalue reference, so can
    // be called with a temporary,
    // i.e. addInactiveClause(vector<Lit>{x, ~y, z}). Returns the
    // added clause, but it is owned by the solver and can be removed
    // at any time. The only safe thing to do with this pointer is to
    // give it back to the solver as the reason for failure/pruning
    template <typename veclit> Clause *addInactiveClause(veclit &&c);

    // Variable mode:
    //
    void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

    // Read state:
    //
    lbool   value      (Var x) const;       // The current value of a variable.
    lbool   value      (Lit p) const;       // The current value of a literal.
    lbool   modelValue (Var x) const;       // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
    std::pair<int,int>  cspModelRange(cspvar x) const; // Range in last model
    int     cspModelValue(cspvar x) const;  // Assigned value in last model
    std::pair< std::set<int>, std::set<int> > const&
            cspSetModel(setvar x) const;    // lb,ub in last model

    int     nAssigns   ()      const;       // The current number of assigned literals.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nVars      ()      const;       // The current number of variables.
    int     nCSPVars   ()      const;       // The number of CSP variables
    int     nConstraints()     const;       // The number of original constraints
    int     nLiveVars()        const;       // Number of vars not fixed at the root

    // Only useful during solving
    int      decisionLevel    ()      const; // Gives the current decisionlevel.
    Lit      decisionAtLevel(int lvl) const; // current decision at given level

    int      varLevel  (Var x) const; // the level where x was assigned. Assumes (but does not assert) value(x) != l_Undef
    int      varLevel  (Lit l) const; // shortcut for varLevel(var(l))

    // the clause that forced x. This is not const because it may
    // require asking an explainer to generate an explicit clause
    Clause *varReason(Var x);
    // shortcut for varReason(var(l)). the clause may have forced ~l
    // instead of l
    Clause *varReason(Lit l);

    // user branching heuristics
    std::function<void(std::vector<Lit>&)> user_brancher; // if user sets varbranch == VAR_USER, this must be non-empty and generate a set of candidates
    lbool currentVarPhase(Var x) const; // give out info to user branchers

    // a fairly dangerous API: get the VSIDS heap with a non-const
    // reference. Useful for implementing variants of VSIDS, but the
    // caller has to maintain the invariant that all unassigned
    // decisions variables are in the heap and that it is a heap
    Heap<VarOrderLt>& vsids_heap();

    // more mild version of the above: just get the current vsids
    // score of a variable.
    double var_activity(Var x);

    // user callback to be notified of every learned clause. gets the
    // clause and the backtrack level. The user can modify the clause,
    // as long as it is correct
    // XXX: the great vec<> vs std::vector<> divide. ugh.
    enum clause_callback_result_t { CCB_OK, CCB_INTERRUPT, CCB_MODIFIED };
    using clause_callback_t =
        std::function<clause_callback_result_t(vec<Lit> &, int)>;

    void use_clause_callback(clause_callback_t cb)
    {
        clause_callbacks.push_back(cb);
    }

    using decision_callback_t = std::function<void()>;

    void use_decision_callback(decision_callback_t cb)
    {
        decision_callbacks.push_back(cb);
    }

    /* User callback to replace the part of UP which finds a new
       watch. It returns UP_result which will have either newwatchidx
       >= 1 (an index into the clause) or < 0, indicating the clause
       is failed. It is not allowed to be 0 or 1

       The callback is allowed to propagate. */
    struct UP_result {
        int newwatchidx;
    };
    using UP_callback_t = std::function<UP_result(Lit, const Clause&)>;

    void use_UP_callback(UP_callback_t cb)
    {
        UP_callback = cb;
    }

    // Extra results: (read-only member variable)
    //
    vec<lbool> model;                   // If problem is satisfiable, this vector contains the model (if any).
    vec<Lit>   conflict;                // If problem is unsatisfiable (possibly under assumptions),
                                        // this vector represent the final conflict clause expressed in the assumptions.
    vec<std::pair<int, int> > cspmodel; // for a satisfiable problem, the lb and ub of every csp variable. easier to
                                        // access than through model, but strictly speaking redundant
    vec<std::pair< std::set<int>, std::set<int> > >
                cspsetmodel;            // for a satisfiable problem, the lb and ub of every set variable

    // Mode of operation:
    //
    bool      trace;              // if true, will output all prunings
    bool      debugclauses;       // if true, will debug clauses generates by constraints *in debug mode only*
    bool      learning;           // can be set to false to disable learning
    bool      restarting;         // can be set to false to disable restarting
    double    var_decay;          // Inverse of the variable activity decay factor.                                            (default 1 / 0.95)
    double    clause_decay;       // Inverse of the clause activity decay factor.                                              (1 / 0.999)
    double    random_var_freq;    // The frequency with which the decision heuristic tries to choose a random variable.        (default 0.02)
    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)
    bool      expensive_ccmin;    // Controls conflict clause minimization.                                                    (default TRUE)
    int       polarity_mode;      // Controls which polarity the decision heuristic chooses. See enum below for allowed modes. (default polarity_false)
    int       verbosity;          // Verbosity level. 0=silent, 1=some progress report                                         (default 0)
    bool      phase_saving;          // standard phase saving: remember last seen polarity in the trail and use it for decisions
    bool      solution_phase_saving; // solution phase saving (for opt problems): remember polarity in last solution and use that
    bool      allow_clause_dbg;   // set to 0 when the solver is cloned to avoid infinite recursion

    bool      interrupt_requested{false}; // true if a callback asked us to stop
    int64_t   conflict_lim{-1};           // stop after this many conflicts

    BranchHeuristic varbranch;
    ValBranchHeuristic valbranch;

    enum { polarity_true = 0, polarity_false = 1, polarity_user = 2, polarity_rnd = 3 };

    // Statistics: (read-only member variable)
    //
    uint64_t starts, decisions, rnd_decisions, propagations, conflicts;
    uint64_t clauses_literals, learnts_literals, max_literals, tot_literals;

public:
    // Interface to propagators
    //
    // Note below that both clauses and explainers are owned by the
    // constraint. It can hand over responsibility for clauses to the
    // solver via addInactiveClause(), but not for explainers
    void     debugclause      (Clause *from, cons *c);                                 // check whether clause *from causes a failure given constraint c
    void     uncheckedEnqueue (Lit p, Clause* from = nullptr);                         // Enqueue a literal. Assumes value of literal is undefined.
    void     uncheckedEnqueueDeferred (Lit p, explainer* from);                        // as above, but with a deferred explanation
    bool     enqueue          (Lit p, Clause* from = nullptr);                         // Test if fact 'p' contradicts current state, enqueue otherwise.
    bool     enqueueDeferred  (Lit p, explainer* from = nullptr);                      // As above, but with a deferred explanation
    Clause*  enqueueFill      (Lit p, vec<Lit>& ps);                                   // Create an inactive clause that includes p and ps and force p

    // A somewhat brittle interface: a non-monotone propagator can
    // discover during pruning that it can in fact backprune a literal
    // (i.e., prune it to a decision level earlier than the current
    // one). It can signal this using the following function. The
    // propagator then has to immediately return the clause 'from' to
    // the propagator, as if it is a conflict. It is not safe to do
    // anything else after calling this.
    //
    // There is no 'deferred' version of this, as the solver would
    // turn around and immediately ask for an explicit explanation
    // anyway.
    //
    // XXX: An alternate, safer interface: have the propagator return
    // a variant<OK, Fail(Clause), Backprune(Lit, Clause)> or
    // something to that effect. But that would require changing all
    // propagators (a fairly mechanical but extensive change).
    void     nonMonotoneEnqueue(Lit p, Clause *from);

    Clause*  propagate        ();                                                      // Perform unit and constraint propagation.
                                                                                       // Returns conflicting clause or NULL

    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     cancelUntil      (int level);                                             // Backtrack until a certain level.

    int cspvarmax(cspvar x);              // get the current max of x
    int cspvarmin(cspvar x);              // get the current min of x
    int cspvardsize(cspvar x);            // get the current domain size of x (which may be smaller than max - min + 1)
    int cspvaromax(cspvar x);             // get the initial max of x
    int cspvaromin(cspvar x);             // get the initial min of x
    Var cspvareqi(cspvar x, int d);       // get the propositional var representing x = d
    Var cspvareqiunsafe(cspvar x, int d); // get the propositional var representing x = d, for the brave
    Var cspvarleqi(cspvar x, int d);      // get the propositional var representing x <= d
    Var cspvarleqiunsafe(cspvar x, int d);// get the propositional var representing x <= d, for the brave
    int setvarumin(setvar x) const;       // get the smallest element in the universe of x
    int setvarumax(setvar x) const;       // get the smallest element in the universe of x
    Var setvarini(setvar x, int d);       // get the propositional var representing d in x
    cspvar setvarcard(setvar x);          // get the cspvar representing the cardinality of x

    // functions that need to be tested, so cannot be private
    void reduce_var(vec<Lit>& ps);        // ensure the domains of the cspvars
                                          // in ps are described minimally
protected:

    // Helper structures:
    //
    friend class VarFilter;
    struct VarFilter {
        const Solver& s;
        VarFilter(const Solver& _s) : s(_s) {}
        bool operator()(Var v) const { return toLbool(s.assigns[v]) == l_Undef && s.decision_var[v]; }
    };

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    vec<Clause*>        clauses;          // List of problem clauses.
    vec<Clause*>        learnts;          // List of learnt clauses.
    vec<Clause*>        inactive;         // List of non-propagating clauses.
    vec<cons*>          conses;           // List of problem constraints.
    vec<consqueue>      consqs;           // Stubs for constraints in the propagation queue
    consqueue*          prop_queue;       // propagation queue (aggregate for all priorities)
    double              cla_inc;          // Amount to bump next clause with.
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.

    struct watch {
      Clause *c{nullptr};
      Lit block{lit_Undef};

      bool operator==(const watch &w) const { return c == w.c; }
      bool operator!=(const watch &w) const { return c != w.c; }
    };

    vec<vec<watch>>     watches;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    vec<vec<watch>>     binwatches;       // same, for binary clauses
    vec<vec<wake_stub>> wakes_on_lit;     // 'wakes_on_lit[var(lit)]' is a list of csp constraints that wake when var is set
    vec<vec<int>>       sched_on_lit;     // 'wakes_on_lit[var(lit)]' is a list of csp constraints that wake when var is set

    vec<char>           assigns;          // The current assignments (lbool:s stored as char:s).
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision_var;     // Declares if a variable is eligible for selection in the decision heuristic.
    vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
    vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<explanation_ptr> reason;          // 'reason[var]' is the clause that implied the variables current value, or 'NULL' if none.
    vec<int>            level;            // 'level[var]' contains the level at which the assignment was made.
    vec<lbool>          phase;            // backjumped over as this polarity
    int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double              random_seed;      // Used by the random variable selection.
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.

    // csp stuff
    vec<cspvar_fixed>   cspvars;             // the fixed data for each cspvar
    vec<domevent>       events;              // the csp event that a literal corresponds to
    vec<setvar_data>    setvars;             // the fixed data for each setvar
    vec<setevent>       setevents;           // the set csp event that a literal corresponds to
    size_t              backtrackable_size;  // How much we need to copy
    size_t              backtrackable_cap;   // How much backtrackable memory is allocated
    vec<void*>          backtrackable_space; // per-level copies of backtrackable data
    void*               current_space;       // All backtrackable data are pointers into this

    cons*               active_constraint;   // the constraint currently propagating.

    std::vector<clause_callback_t> clause_callbacks; // all clause callbacks
    std::vector<decision_callback_t> decision_callbacks; // all decision callbacks
    UP_callback_t       UP_callback;

    // names. for tracing output etc
    std::vector<std::string>  varnames;
    std::vector<std::string>  cspvarnames;
    std::vector<std::string>  setvarnames;

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            analyze_explbuffer;
    vec<Lit>            add_tmp;

    std::vector<char> reduce_var_seen;
    std::vector<int> reduce_var_min;
    std::vector<int> reduce_var_max;
    std::vector<char> reduce_var_asgn;
    std::vector<cspvar> reduce_var_toclear;

    std::vector<Lit> user_candidates;

    // Main internal methods:
    //
    void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
    Lit      pickBranchLit     (int polarity_mode, double random_var_freq);            // Return the next decision variable.
    Lit      pickBranchLitVSIDS(int polarity_mode, double random_var_freq);            // ...using VSIDS
    Lit      pickBranchLitLex  ();                                                     // ...Lex order of the CSP vars
    Lit      pickBranchLitDom  ();                                                     // ...mindom among csp vars
    Lit      pickBranchLitFrom (cspvar x);                                             // Branch within a chosen variable
    void     analyze          (Clause* confl, vec<Lit>& out_learnt, int& out_btlevel); // (bt = backtrack)
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')
    lbool    search           (int nof_conflicts, double * nof_learnts);               // Search for a given number of conflicts.
    void     reduceDB         ();                                                      // Reduce the set of learnt clauses.
    void     gcInactive       ();                                                      // Collect inactive and unlocked clauses
    void     removeSatisfied  (vec<Clause*>& cs);                                      // Shrink 'cs' to contain only non-satisfied clauses.
    void     uncheckedEnqueue_common(Lit p, explanation_ptr from);                     // common stuff done by both uE(Lit, Clause*) and uE(Lit, explainer*)
    void     uncheckedEnqueue_np(Lit p, explanation_ptr from);                         // uncheckedEnqueue with no CSP propagation
    Clause*  propagate_inner  ();                                                      // Perform unit propagation, wake propagators,
                                                                                       // schedule propagators. Returns conflicting clause or NULL
    Clause*  propagate_wakes  (Lit p, const vec<wake_stub>& wakes);                    // wake all constraints that are registered in this list
    Clause*  explicit_reason  (Lit p);                                                 // If a literal has a Clause as reason, return that. Otherwise, use the explainer to create
                                                                                       // a clause, change the reason to that new clause and return that.
    bool     withinBudget     () const;                                                // if true, we have exceeded the budget

    // Constraint queue
    //
    void     ensure_can_schedule(cons *c);             // Create a consqueue for c
    void     schedule(int idx);                        // Put consqs[idx].c in its queue
    int      first_scheduled();                        // Head of the queue, <0 if queue is empty
    void     unschedule(int idx);                      // Remove consqs[idx].c from the queue
    void     reset_queue();                            // Remove everything from the queue

    // Maintaining Variable/Clause activity:
    //
    void     varDecayActivity ();                      // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void     varBumpActivity  (Var v);                 // Increase a variable with the current 'bump' value.
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.

    // Operations on clauses:
    //
    void     attachClause     (Clause& c);             // Attach a clause to watcher lists.
    void     detachClause     (Clause& c);             // Detach a clause to watcher lists.
    void     removeClause     (Clause& c);             // Detach and free a clause.
    bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    // Misc:
    //
    uint32_t abstractLevel    (Var x) const; // Used to represent an abstraction of sets of decision levels.
    double   progressEstimate ()      const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    // Debug:
    void     printLit         (Lit l);
    template<class C>
    void     printClause      (const C& c);
    void     verifyModel      ();
    void     checkLiteralCount();

    friend std::ostream& operator<<(std::ostream&, lit_printer);

public:
    std::vector<int> debug_solution;
    std::vector<lbool> debug_solution_lits;
    void check_debug_solution(Lit p, explanation_ptr from);

    // Static helpers:
    //

    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size); }

    void    setrandomseed(double s) { random_seed = s; }
};

//=================================================================================================
// Implementation of inline methods:

inline
void* Solver::get(btptr p)
{
  return ((char*)current_space)+p.offset;
}

template <typename T>
inline
T& Solver::deref(btptr p)
{
  return *(T*)get(p);
}

template <typename T>
inline
T* Solver::deref_array(btptr p)
{
  return (T*)get(p);
}

inline
domevent Solver::event(Lit p) const
{
  if( p != lit_Undef ) return events[toInt(p)];
  else return domevent();
}

inline
setevent Solver::sevent(Lit p) const
{
  if( p != lit_Undef ) return setevents[toInt(p)];
  else return setevent();
}

inline void Solver::insertVarOrder(Var x) {
    if (!order_heap.inHeap(x) && decision_var[x]) order_heap.insert(x); }

inline void Solver::varDecayActivity() { var_inc *= var_decay; }
inline void Solver::varBumpActivity(Var v) {
    if ( (activity[v] += var_inc) > 1e100 ) {
        // Rescale:
        for (int i = 0; i < nVars(); i++)
            activity[i] *= 1e-100;
        var_inc *= 1e-100; }

    // Update order_heap with respect to new activity:
    if (order_heap.inHeap(v))
        order_heap.decrease(v); }

inline void Solver::claDecayActivity() { cla_inc *= clause_decay; }
inline void Solver::claBumpActivity (Clause& c) {
        if ( (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++)
                learnts[i]->activity() *= 1e-20;
            cla_inc *= 1e-20; } }

inline Heap<Solver::VarOrderLt>& Solver::vsids_heap() { return order_heap; }

inline double Solver::var_activity(Var x) { return activity[x]; }

inline bool Solver::enqueue(Lit p, Clause* from)
{
    return value(p) != l_Undef ? value(p) != l_False
                               : (uncheckedEnqueue(p, from), true);
}

inline bool Solver::enqueueDeferred(Lit p, explainer* from)
{
    if (value(p) != l_Undef)
        return value(p) != l_False;
    uncheckedEnqueueDeferred(p, from);
    return true;
}

inline bool Solver::locked(const Clause &c) const {
    auto c0r = reason[var(c[0])];
    if (!c0r.has<Clause>())
        return false;
    auto *cl = c0r.get<Clause>();
    return cl == &c && value(c[0]) == l_True;
}

inline void     Solver::newDecisionLevel()  {
  trail_lim.push(trail.size());
  int dlvl = decisionLevel();
  while( backtrackable_space.size() < decisionLevel() )
    backtrackable_space.push(0L);
  if( !backtrackable_space[dlvl-1] ) {
    backtrackable_space[dlvl-1] = malloc(backtrackable_size);
    if(! backtrackable_space[dlvl-1] )
      throw std::bad_alloc();
  }
  std::memcpy(backtrackable_space[dlvl-1], current_space,
              backtrackable_size);
}

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }

inline Lit      Solver::decisionAtLevel(int lvl) const
{
    assert(lvl >= 1 && lvl <= decisionLevel());
    auto pos = trail_lim[lvl - 1];
    if (lvl < decisionLevel() && pos == trail_lim[lvl])
        return lit_Undef;
    if (lvl == decisionLevel() && pos == trail.size())
        return lit_Undef;
    return trail[pos];
}

inline int Solver::varLevel(Var x) const { return level[x]; }
inline int Solver::varLevel(Lit l) const { return level[var(l)]; }

inline Clause *Solver::varReason(Var x) { return explicit_reason(Lit(x)); }
inline Clause *Solver::varReason(Lit l) { return explicit_reason(l); }

inline uint32_t Solver::abstractLevel (Var x) const   { return 1 << (level[x] & 31); }
inline lbool    Solver::value         (Var x) const   { return toLbool(assigns[x]); }
inline lbool    Solver::value         (Lit p) const   { return toLbool(assigns[var(p)]) ^ sign(p); }
inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return assigns.size(); }
inline int      Solver::nCSPVars      ()      const   { return cspvars.size(); }
inline int      Solver::nConstraints  ()      const   { return conses.size(); }
inline int Solver::nLiveVars() const
{
    if (decisionLevel() > 0)
        return nVars() - trail_lim[0];
    else
        return nVars() - trail.size();
}
inline void     Solver::setPolarity   (Var v, bool b) { polarity    [v] = (char)b; }
inline void     Solver::setDecisionVar(Var v, bool b) { decision_var[v] = (char)b; if (b) { insertVarOrder(v); } }
inline bool     Solver::solve         ()              { vec<Lit> tmp; return solve(tmp); }
inline lbool    Solver::solveBudget   ()              { vec<Lit> tmp; return solveBudget(tmp); }
inline bool     Solver::okay          ()      const   { return ok; }
inline bool     Solver::withinBudget  ()      const   { return conflict_lim < 0 || conflicts < static_cast<uint64_t>(conflict_lim); }

template <typename veclit>
inline Clause *Solver::addInactiveClause(veclit &&ps) {
  for (Lit l : ps)
    assert(var(l) != var_Undef);
  Clause *r = Clause_new(ps);
  addInactiveClause(r);
  return r;
}

template <typename veclit> void Solver::addClause(veclit &&v)
{
    add_tmp.clear();
    using std::begin;
    using std::end;
    std::copy(begin(v), end(v), std::back_inserter(add_tmp));
    addClause(add_tmp);
}

//=================================================================================================
// Debug + etc:

#ifdef _MSC_VER
#define reportf(format, ...) ( fflush(stdout), fprintf(stdout, format, ## __VA_ARGS__), fflush(stdout) )
#else
#define reportf(format, args...) ( fflush(stdout), fprintf(stdout, format, ## args), fflush(stdout) )
#endif

static inline void logLit(FILE* f, Lit l)
{
    fprintf(f, "%sx%d", sign(l) ? "~" : "", var(l)+1);
}

static inline void logLits(FILE* f, const vec<Lit>& ls)
{
    fprintf(f, "[ ");
    if (ls.size() > 0){
        logLit(f, ls[0]);
        for (int i = 1; i < ls.size(); i++){
            fprintf(f, ", ");
            logLit(f, ls[i]);
        }
    }
    fprintf(f, "] ");
}

static inline const char* showBool(bool b) { return b ? "true" : "false"; }


// Just like 'assert()' but expression will be evaluated in the release version as well.
static inline void check(bool expr) { assert(expr); }


inline void Solver::printLit(Lit l)
{
    reportf("%s%d:%c", sign(l) ? "-" : "", var(l)+1, value(l) == l_True ? '1' : (value(l) == l_False ? '0' : 'X'));
}


template<class C>
inline void Solver::printClause(const C& c)
{
    for (int i = 0; i < c.size(); i++){
        printLit(c[i]);
        fprintf(stderr, " ");
    }
}

inline
std::ostream& operator<<(std::ostream& os, var_printer vp)
{
  std::string const& vn = vp._s.getVarName(vp._v);
  if( vn.empty() )
    os << vp._v+1;
  else
    os << vn;
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, cspvar_printer vp)
{
  std::string const& vn = vp._s.getCSPVarName(vp._v);
  if( vn.empty() )
    os << "x" << vp._v.id();
  else
    os << vn;
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, setvar_printer vp)
{
  std::string const& vn = vp._s.getSetVarName(vp._v);
  if( vn.empty() )
    os << "s" << vp._v.id();
  else
    os << vn;
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, lit_printer lp)
{
  Solver &s(lp._s);
  Lit l = lp._p;
  if (var(l) == var_Undef) {
    os << "lit_Undef";
    return os;
  }
  domevent pe = s.event(l);
  if( noevent(pe) )
    os << (sign(l) ? "-" : "") << var_printer(s, var(l));
  else
    os << domevent_printer(s, pe);
  os << ":"
     << (s.value(l) == l_True ? '1' : (s.value(l) == l_False ? '0' : 'X'));
  if( s.value(l) != l_Undef )
    os << '(' << s.level[var(l)] << ')';
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, domain_as_range d)
{
  cspvar x = d._x;
  Solver &s = d._s;
  os << "[" << x.min(s) << ".." << x.max(s) << "]";
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, domain_as_set d)
{
  cspvar x = d._x;
  Solver &s = d._s;
  os << "{" << x.min(s);
  int beg = x.min(s), end = beg;
  for(int i = x.min(s)+1; i <= x.max(s); ++i) {
    if( !x.indomain(s, i) ) {
      if( end < i - 1 ) continue;
      if( end > beg )
        os << ".." << end;
      os << ",";
      beg = i;
    } else {
      if( beg > end ) {
        os << i;
        beg = i;
      }
      end = i;
    }
  }
  if( end > beg )
    os << ".." << end;
  os << "}";
  return os;
}

inline std::pair<int, int> Solver::cspModelRange(cspvar x) const
{
  return cspmodel[x._id];
}

inline int Solver::cspModelValue(cspvar x) const
{
  int id = x._id;
  if( cspmodel[id].first != cspmodel[id].second )
    throw unassigned_var(*this, x);
  return cspmodel[id].first;
}

inline
std::pair< std::set<int>, std::set<int> > const&
Solver::cspSetModel(setvar x) const
{
  return cspsetmodel[x._id];
}


//==================================================
// cspvar inlines

inline int Solver::cspvarmax(cspvar x)
{
  return cspvars[x._id].max;
}

inline int Solver::cspvarmin(cspvar x)
{
  return cspvars[x._id].min;
}

inline int Solver::cspvardsize(cspvar x)
{
  return cspvars[x._id].dsize;
}

inline int Solver::cspvaromin(cspvar x)
{
  return cspvars[x._id].omin;
}

inline int Solver::cspvaromax(cspvar x)
{
  return cspvars[x._id].omax;
}

inline Var Solver::cspvareqi(cspvar x, int d)
{
  return cspvars[x._id].eqi(d);
}

inline Var Solver::cspvareqiunsafe(cspvar x, int d)
{
  return cspvars[x._id].eqiUnsafe(d);
}

inline Var Solver::cspvarleqi(cspvar x, int d)
{
  return cspvars[x._id].leqi(d);
}

inline Var Solver::cspvarleqiunsafe(cspvar x, int d)
{
  return cspvars[x._id].leqiUnsafe(d);
}

inline Var Solver::setvarumin(setvar x) const
{
  return setvars[x._id].min;
}

inline Var Solver::setvarumax(setvar x) const
{
  return setvars[x._id].max;
}

inline Var Solver::setvarini(setvar x, int d)
{
  return setvars[x._id].ini(d);
}

inline cspvar Solver::setvarcard(setvar x)
{
  return setvars[x._id]._card;
}

inline bool cspvar::indomain(Solver &s, int d) const
{
  Var xd = eqi(s, d);
  return xd != var_Undef && s.value( xd ) != l_False;
}

inline bool cspvar::indomainUnsafe(Solver &s, int d) const
{
  Var xd = eqiUnsafe(s, d);
  return s.value( xd ) != l_False;
}

inline int cspvar::min(Solver &s) const
{
  return s.cspvarmin(*this);
}

inline int cspvar::max(Solver &s) const
{
  return s.cspvarmax(*this);
}

inline int cspvar::domsize(Solver &s) const
{
  return s.cspvardsize(*this);
}

inline int cspvar::omin(Solver &s) const
{
  return s.cspvaromin(*this);
}

inline int cspvar::omax(Solver &s) const
{
  return s.cspvaromax(*this);
}

inline Var cspvar::eqi(Solver &s, int d) const
{
  return s.cspvareqi(*this, d);
}

inline Var cspvar::eqiUnsafe(Solver &s, int d) const
{
  return s.cspvareqiunsafe(*this, d);
}

inline Var cspvar::leqi(Solver &s, int d) const
{
  return s.cspvarleqi(*this, d);
}

inline Var cspvar::leqiUnsafe(Solver &s, int d) const
{
  return s.cspvarleqiunsafe(*this, d);
}

inline Lit cspvar::r_geq(Solver &s, int d) const
{
  return Lit( leqi(s, d-1) );
}

inline Lit cspvar::r_leq(Solver &s, int d) const
{
  return ~Lit( leqi(s, d) );
}

inline Lit cspvar::r_eq(Solver &s, int d) const
{
  return ~Lit( eqi(s, d) );
}

inline Lit cspvar::r_eq(Solver &s) const
{
  Lit l = ~Lit( eqi(s, min(s)));
  assert( s.value(l) == l_False );
  return l;
}

inline Lit cspvar::r_neq(Solver &s, int d) const
{
  return Lit( eqi(s, d) );
}

inline Lit cspvar::r_min(Solver &s) const
{
  return r_geq(s, min(s));
}

inline Lit cspvar::r_max(Solver &s) const
{
  return r_leq(s, max(s));
}

inline Lit cspvar::e_geq(Solver &s, int d) const
{
  return ~Lit( leqi(s, d-1) );
}

inline Lit cspvar::e_leq(Solver &s, int d) const
{
  return Lit( leqi(s, d) );
}

inline Lit cspvar::e_eq(Solver &s, int d) const
{
  return Lit( eqi(s, d) );
}

inline Lit cspvar::e_neq(Solver &s, int d) const
{
  return ~Lit( eqi(s, d) );
}

inline Clause *cspvar::remove(Solver &s, int d, Clause *c)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  if( s.value(xd) == l_True ) return throw_if_null(c);
  s.uncheckedEnqueue( ~Lit(xd), c);
  return 0L;
}

template<typename V>
inline Clause *cspvar::remove(Solver &s, int d, V& ps)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  Clause *r = Clause_new(ps, false, ~Lit(xd));
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r);
  return 0L;
}

template<typename V>
inline Clause *cspvar::removef(Solver &s, int d, V& ps)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  ps.push( ~Lit(xd) );
  Clause *r = Clause_new(ps, false, ~Lit(xd));
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r);
  return 0L;
}

inline Clause *cspvar::setmin(Solver &s, int d, Clause *c)
{
  Var xd = leqi(s, d-1);
  if( xd == var_Undef ) {
    if( d <= max(s) ) return 0L;
    return throw_if_null(c);
  }
  if( s.value(xd) == l_False ) return 0L;
  if( s.value(xd) == l_True ) return throw_if_null(c);
  s.uncheckedEnqueue( ~Lit(xd), c);
  return 0L;
}

template<typename V>
inline Clause *cspvar::setmin(Solver &s, int d, V& ps)
{
  Var xd = leqi(s, d-1);
  if( xd == var_Undef ) {
    if( d <= max(s) ) return 0L;
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_False ) return 0L;
  Clause *r = Clause_new(ps, false, ~Lit(xd) );
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r);
  return 0L;
}

template<typename V>
inline Clause *cspvar::setminf(Solver &s, int d, V& ps)
{
  Var xd = leqi(s, d-1);
  if( xd == var_Undef ) {
    if( d <= max(s) ) return 0L;
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_False ) return 0L;
  ps.push( ~Lit(xd) );
  Clause *r = Clause_new(ps, false, ~Lit(xd) );
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r);
  return 0L;
}

inline Clause *cspvar::setmax(Solver &s, int d, Clause *c)
{
  Var xd = leqi(s, d);
  if( xd == var_Undef ) {
    if( d >= min(s) ) return 0L;
    return throw_if_null(c);
  }
  if( s.value(xd) == l_True ) return 0L;
  if( s.value(xd) == l_False ) return throw_if_null(c);
  s.uncheckedEnqueue( Lit(xd), c);
  return 0L;
}

template<typename V>
inline Clause *cspvar::setmax(Solver &s, int d, V& ps)
{
  Var xd = leqi(s, d);
  if( xd == var_Undef ) {
    if( d >= min(s) ) return 0L;
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  Clause *r = Clause_new(ps, false, Lit(xd) );
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r);
  return 0L;
}

template<typename V>
inline Clause *cspvar::setmaxf(Solver &s, int d, V& ps)
{
  Var xd = leqi(s, d);
  if( xd == var_Undef ) {
    if( d >= min(s) ) return 0L;
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  ps.push( Lit(xd) );
  Clause *r = Clause_new(ps, false, Lit(xd) );
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r);
  return 0L;
}

inline Clause *cspvar::assign(Solver &s, int d, Clause *c)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef )
    return throw_if_null(c);
  if( s.value(xd) == l_True ) return 0L;
  if( s.value(xd) == l_False ) return c;
  s.uncheckedEnqueue( Lit(xd), c );
  return 0L;
}

template<typename V>
inline Clause *cspvar::assign(Solver &s, int d, V& ps)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef ) {
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  Clause *r = Clause_new(ps, false, Lit(xd) );
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r );
  return 0L;
}

template<typename V>
inline Clause *cspvar::assignf(Solver &s, int d, V& ps)
{
  Var xd = eqi(s, d);
  if( xd == var_Undef ) {
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  ps.push( Lit(xd) );
  Clause *r = Clause_new(ps, false, Lit(xd) );
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r );
  return 0L;
}

inline bool operator==(cspvar x1, cspvar x2)
{
  return x1.id() == x2.id();
}

inline bool operator<(cspvar x1, cspvar x2)
{
  return x1.id() < x2.id();
}

inline
int setvar::umin(Solver &s) const {
  return s.setvarumin(*this);
}

inline
int setvar::umax(Solver &s) const {
  return s.setvarumax(*this);
}

inline
cspvar setvar::card(Solver &s) const {
  return s.setvarcard(*this);
}

inline
int setvar::ini(Solver& s, int d) const {
  return s.setvarini(*this, d);
}

inline
bool setvar::includes(Solver &s, int d) const {
  Var xd = ini(s, d);
  if( xd == var_Undef ) return false;
  else return (s.value(xd) == l_True);
}

inline
bool setvar::excludes(Solver &s, int d) const {
  Var xd = ini(s, d);
  if( xd == var_Undef ) return true;
  else return (s.value(xd) == l_False);
}

inline Clause *setvar::exclude(Solver &s, int d, Clause *c)
{
  Var xd = ini(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  if( s.value(xd) == l_True ) return c;
  s.uncheckedEnqueue( ~Lit(xd), c );
  return 0L;
}

inline Clause *setvar::exclude(Solver &s, int d, vec<Lit>& ps)
{
  Var xd = ini(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  Clause *r = Clause_new(ps, false, ~Lit(xd) );
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r );
  return 0L;
}

inline Clause *setvar::excludef(Solver &s, int d, vec<Lit>& ps)
{
  Var xd = ini(s, d);
  if( xd == var_Undef ) return 0L;
  if( s.value(xd) == l_False ) return 0L;
  ps.push( ~Lit(xd) );
  Clause *r = Clause_new(ps, false, ~Lit(xd) );
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_True ) return r;
  s.uncheckedEnqueue( ~Lit(xd), r );
  return 0L;
}

inline Clause *setvar::include(Solver &s, int d, Clause *c)
{
  Var xd = ini(s, d);
  if( xd == var_Undef )
    throw_if_null(c);
  if( s.value(xd) == l_True ) return 0L;
  if( s.value(xd) == l_False ) return c;
  s.uncheckedEnqueue( Lit(xd), c );
  return 0L;
}

inline Clause *setvar::include(Solver &s, int d, vec<Lit>& ps)
{
  Var xd = ini(s, d);
  if( xd == var_Undef ) {
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  Clause *r = Clause_new(ps, false, Lit(xd) );
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r );
  return 0L;
}

inline Clause *setvar::includef(Solver &s, int d, vec<Lit>& ps)
{
  Var xd = ini(s, d);
  if( xd == var_Undef ) {
    Clause *r = Clause_new(ps);
    s.addInactiveClause(r);
    return r;
  }
  if( s.value(xd) == l_True ) return 0L;
  ps.push( Lit(xd) );
  Clause *r = Clause_new(ps, false, Lit(xd) );
  ps.pop();
  s.addInactiveClause(r);
  if( s.value(xd) == l_False ) return r;
  s.uncheckedEnqueue( Lit(xd), r );
  return 0L;
}

inline bool operator==(setvar x1, setvar x2)
{
  return x1.id() == x2.id();
}

//==================================================
//

// push p to ps iff var(p) >= 0. Return 1 if we did push, 0 otherwise
template<typename V>
inline int pushifdef(V& ps, Lit p)
{
  if( var(p) >= 0 ) {
    ps.push(p);
    return 1;
  }
  return 0;
}

struct push_temp_p {
  vec<Lit> & _ps;
  Lit _p;
  push_temp_p(vec<Lit>& ps, Lit p) : _ps(ps), _p(p) {
    if( _p != lit_Undef ) _ps.push(_p);
  }
  ~push_temp_p() { if(_p != lit_Undef ) _ps.pop(); }
};

#ifdef _MSC_VER
#define PUSH_TEMP(x, y) push_temp_p BOOST_PP_CAT(pp, __LINE__)(x, y)
#else
#define PUSH_TEMP(x, y) push_temp_p BOOST_PP_CAT(pp, __LINE__) \
  __attribute__((unused)) (x,y)
#endif

//==================================================
// more usable backtrackables
template <typename T> struct backtrackable {
  static_assert(std::is_pod_v<T>, "Only POD types can be backtrackable");
  Solver &s;
  btptr p;

  backtrackable(const backtrackable &) = default;
  backtrackable(Solver &s, const T &t = T{})
      : s(s), p(s.alloc_backtrackable(sizeof(T))) {
    T &st = s.deref<T>(p);
    new (&st) T(t);
  }

  T &operator*() { return s.deref<T>(p); }
  const T &operator*() const { return s.deref<T>(p); }
  T *operator->() { return &s.deref<T>(p); }
  const T *operator->() const { return &s.deref<T>(p); }
};

//==================================================
// a trick to avoid branching
// returns a1 if w >= 0, otherwise a2

inline
int select(int w, int a1, int a2)
{
  int mask = w >> (sizeof(int)*8-1);
  return - (~mask*a1) - (mask*a2);
}

//==================================================
// exception stuff

inline
const char *unassigned_var::what() const throw()
{
  const char s[] = "expected var x%d to be assigned";
  static const int l = strlen(s);
  static char exc[2*sizeof(s)+1];
  snprintf(exc, 2*l, s, _x.id());
  return exc;
}

inline
const char *undefined_literal::what() const throw()
{
  const char s[] = "literal x%d %s %d is undefined";
  static const int l = strlen(s);
  static char exc[2*sizeof(s)+1];
  snprintf(exc, 2*l, s, _e.x.id(), opstring(_e.type), _e.d);
  return exc;
}

} // namespace minicsp

//=================================================================================================
#endif
