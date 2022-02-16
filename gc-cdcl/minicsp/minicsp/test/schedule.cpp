/*************************************************************************
minicsp

Copyright 2010--2011 George Katsirelos

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

#include <iostream>
#include "test.hpp"
#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"

using namespace std;

namespace {
  typedef map<size_t, vector<domevent> > action_schedule;

  class test_cons : public cons
  {
  public:
    int _id;
    vector<int> * _wakes; // if NULL, just fail on the first invocation

    vector<cspvar> _d;
    vector<cspvar> _l;
    vector<cspvar> _u;
    vector<cspvar> _f;

    action_schedule _actions;

    test_cons(Solver &s,
              int id,
              vector<int>* wakes,
              vector<cspvar> const& d,
              vector<cspvar> const& l,
              vector<cspvar> const& u,
              vector<cspvar> const& f,
              action_schedule const & actions,
              int priority) :
      _id(id), _wakes(wakes),
      _d(d), _l(l), _u(u), _f(f),
      _actions(actions)
    {
      set_priority(priority);
      for(size_t i = 0; i != _d.size(); ++i)
        s.schedule_on_dom(_d[i], this);
      for(size_t i = 0; i != _l.size(); ++i)
        s.schedule_on_lb(_l[i], this);
      for(size_t i = 0; i != _u.size(); ++i)
        s.schedule_on_ub(_u[i], this);
      for(size_t i = 0; i != _f.size(); ++i)
        s.schedule_on_fix(_f[i], this);
    }

    void add_all_reasons(Solver &s, vec<Lit>& ps, vector<cspvar>& v,
                         cspvar except)
    {
      for(size_t i = 0; i != v.size(); ++i) {
        if( v[i] == except ) continue;
        pushifdef( ps, v[i].r_min(s) );
        pushifdef( ps, v[i].r_max(s) );
        for( int q = v[i].min(s)+1; q < v[i].max(s); ++q)
          if( !v[i].indomain(s, q) )
            ps.push( v[i].r_neq(s, q));
      }
    }

    Clause *propagate(Solver& s)
    {
      if(!_wakes) {
        vec<Lit> ps;
        cspvar x;
        add_all_reasons(s, ps, _d, x);
        add_all_reasons(s, ps, _l, x);
        add_all_reasons(s, ps, _u, x);
        add_all_reasons(s, ps, _f, x);
        Clause *r = Clause_new(ps);
        s.addInactiveClause(r);
        return r;
      }

      (*_wakes).push_back(_id);
      action_schedule::const_iterator i =
        _actions.find(_wakes->size());
      if( i == _actions.end() ) return 0L;

      vec<Lit> ps;
      vector<domevent> const& vde = i->second;
      for(size_t q = 0; q != vde.size(); ++q) {
        domevent de = vde[q];
        ps.clear();

        add_all_reasons(s, ps, _d, de.x);
        add_all_reasons(s, ps, _l, de.x);
        add_all_reasons(s, ps, _u, de.x);
        add_all_reasons(s, ps, _f, de.x);

        Lit p;
        switch( de.type ) {
        case domevent::EQ: p = de.x.e_eq(s, de.d); break;
        case domevent::NEQ: p = de.x.e_neq(s, de.d); break;
        case domevent::LEQ: p = de.x.e_leq(s, de.d); break;
        case domevent::GEQ: p = de.x.e_geq(s, de.d); break;
        case domevent::NONE: assert(0);
        }
        DO_OR_RETURN(s.enqueueFill(p, ps));
      }
      return 0L;
    }

    void clone(Solver &other)
    {
      cons *con = new test_cons(other, _id, 0L, _d, _l, _u, _f,
                                _actions, 0);
      other.addConstraint(con);
    }
  };

  void post_test(Solver &s,
                 int id,
                 vector<int>& wakes,
                 vector<cspvar> const& d,
                 vector<cspvar> const& l,
                 vector<cspvar> const& u,
                 vector<cspvar> const& f)
  {
    action_schedule a;
    cons *con = new test_cons(s, id, &wakes, d, l, u, f, a, 0);
    s.addConstraint(con);
  }

  void post_test(Solver &s,
                 int id,
                 vector<int>& wakes,
                 vector<cspvar> const& d,
                 vector<cspvar> const& l,
                 vector<cspvar> const& u,
                 vector<cspvar> const& f,
                 action_schedule const & actions,
                 int priority = 0)
  {
    cons *con = new test_cons(s, id, &wakes, d, l, u, f, actions, priority);
    s.addConstraint(con);
  }

  bool compare_events(vector<int> const &wakes,
                      int *expected)
  {
    size_t i = 0;
    for(; i != wakes.size() && expected[i] >= 0 ; ++i)
      if(expected[i] != wakes[i])
        return false;
    if( expected[i] >= 0 || i < wakes.size() )
      return false;
    return true;
  }

  // a prop wakes
  void schedule01()
  {
    Solver s;
    vector<cspvar> d, l, u, f;
    d = s.newCSPVarArray(2, 5, 10);

    vector<int> wakes;

    post_test(s, 1, wakes, d, l, u, f);

    int exp0[] = { -1 };
    assert( !s.propagate() );
    assert( compare_events(wakes, exp0) );

    int exp1[] = { 1, -1 };
    s.newDecisionLevel();
    d[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp1) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule01);

  // a prop wakes in all event combinations
  void schedule02()
  {
    Solver s;
    vector<cspvar> d, l, u, f;
    d = s.newCSPVarArray(2, 5, 10);
    l = s.newCSPVarArray(2, 5, 10);
    u = s.newCSPVarArray(2, 5, 10);
    f = s.newCSPVarArray(2, 5, 10);

    vector<int> wakes;

    post_test(s, 1, wakes, d, l, u, f);

    int exp[] = { 1, -1 };
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );

    // lbound change triggers dom schedule
    wakes.clear();
    s.newDecisionLevel();
    d[0].setmin(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // ubound change triggers dom schedule
    wakes.clear();
    s.newDecisionLevel();
    d[0].setmax(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // fix triggers dom schedule
    wakes.clear();
    s.newDecisionLevel();
    d[0].assign(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // dom triggers lbnd schedule
    wakes.clear();
    s.newDecisionLevel();
    l[0].remove(s, 5, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // dom triggers ubnd schedule
    wakes.clear();
    s.newDecisionLevel();
    u[0].remove(s, 10, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // lbnd+ubnd triggers fix schedule
    wakes.clear();
    s.newDecisionLevel();
    f[0].setmin(s, 7, NO_REASON);
    f[0].setmax(s, 7, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // fix triggers dom schedule
    wakes.clear();
    s.newDecisionLevel();
    d[0].assign(s, 7, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // fix triggers lbnd schedule
    wakes.clear();
    s.newDecisionLevel();
    l[0].assign(s, 7, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);

    // fix triggers ubnd schedule
    wakes.clear();
    s.newDecisionLevel();
    u[0].assign(s, 7, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule02);

  // schedule 2 propagators
  void schedule03()
  {
    Solver s;
    vector<cspvar> d0, d1, l, u, f;
    vector<cspvar> x = s.newCSPVarArray(3, 5, 10);
    d0.push_back(x[0]);
    d0.push_back(x[1]);
    d1.push_back(x[1]);
    d1.push_back(x[2]);

    vector<int> wakes;
    post_test(s, 1, wakes, d0, l, u, f);
    post_test(s, 2, wakes, d1, l, u, f);

    int exp1[] = { 1, -1 };
    wakes.clear();
    s.newDecisionLevel();
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp1) );
    s.cancelUntil(0);

    int exp2[] = { 2, -1 };
    wakes.clear();
    s.newDecisionLevel();
    x[2].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp2) );
    s.cancelUntil(0);

    int exp3[] = { 1, 2, -1 };
    wakes.clear();
    s.newDecisionLevel();
    x[0].remove(s, 6, NO_REASON);
    x[2].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp3) );
    s.cancelUntil(0);

    int exp4[] = { 2, 1, -1 };
    wakes.clear();
    s.newDecisionLevel();
    x[2].remove(s, 6, NO_REASON);
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp4) );
    s.cancelUntil(0);

    wakes.clear();
    s.newDecisionLevel();
    x[1].remove(s, 6, NO_REASON);
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp3) || compare_events(wakes, exp4) );
    s.cancelUntil(0);

  }
  REGISTER_TEST(schedule03);

  // reschedule 1 propagator, no others present
  void schedule04()
  {
    Solver s;
    vector<cspvar> d, l, u, f;
    d = s.newCSPVarArray(2, 5, 10);

    vector<int> wakes;

    action_schedule a;
    a[1].push_back(domevent(d[0], domevent::NEQ, 7));
    post_test(s, 1, wakes, d, l, u, f, a);
    post_eq(s, d[0], d[1], 0);

    int exp0[] = { -1 };
    int exp1[] = { 1, 1, -1 };
    assert( !s.propagate() );
    assert( compare_events( wakes, exp0 ) );

    s.newDecisionLevel();
    d[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp1) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule04);

  // reschedule propagator
  void schedule05()
  {
    Solver s;
    vector<cspvar> d0, d1, l, u, f;
    vector<cspvar> x = s.newCSPVarArray(5, 5, 10);
    d0.assign(x.begin(), x.begin()+2);
    d1.assign(x.begin()+1, x.begin()+3);

    vector<int> wakes;

    action_schedule a;
    a[1].push_back(domevent(d0[0], domevent::NEQ, 7));
    post_test(s, 1, wakes, d0, l, u, f, a);
    post_test(s, 2, wakes, d1, l, u, f);

    post_eq(s, x[0], x[3], 0);
    post_eq(s, x[3], x[2], 0);
    post_eq(s, x[2], x[4], 0);
    post_eq(s, x[4], x[0], 0);

    int exp0[] = { -1 };
    int exp1[] = { 1, 2, 1, -1 };
    assert( !s.propagate() );
    assert( compare_events( wakes, exp0 ) );

    s.newDecisionLevel();
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp1) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule05);

  // check that queue is cleared correctly: schedule p, fail before p
  // is executed, backtrack and make decisions that do not schedule p
  // again. p must not be in wakes.
  void schedule_clear01()
  {
    Solver s;
    vector<cspvar> d, l, u, f;
    d = s.newCSPVarArray(3, 5, 10);
    Var b = s.newVar();

    vector<int> wakes;
    post_test(s, 1, wakes, d, l, u, f);
    post_eq(s, d[0], d[1], 0);
    post_less_re(s, d[0], d[1], 0, Lit(b));

    int exp[] = { -1 };

    s.newDecisionLevel();
    s.enqueue( Lit(b) );
    assert( s.propagate() );
    s.cancelUntil(0);

    s.newDecisionLevel();
    s.enqueue( ~Lit(b) );
    assert( !s.propagate() );
    s.cancelUntil(0);

    assert( compare_events(wakes, exp) );
  }
  REGISTER_TEST(schedule_clear01);

  // check that queue is cleared correctly, take 2: schedule p, fail
  // before p is executed, backtrack and make decisions that schedule
  // p again. p must be in wakes.
  void schedule_clear02()
  {
    Solver s;
    vector<cspvar> d, l, u, f;
    d = s.newCSPVarArray(3, 5, 10);
    Var b = s.newVar();

    vector<int> wakes;
    post_test(s, 1, wakes, d, l, u, f);
    post_eq(s, d[0], d[1], 0);
    post_less_re(s, d[0], d[1], 0, Lit(b));

    int exp[] = { 1, -1 };

    s.newDecisionLevel();
    s.enqueue( Lit(b) );
    assert( s.propagate() );
    s.cancelUntil(0);

    s.newDecisionLevel();
    d[0].remove(s, 7, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);

    assert( compare_events(wakes, exp) );
  }
  REGISTER_TEST(schedule_clear02);


  //--------------------------------------------------
  // multiple priorities

  // 2 priority 0, 2 priority 1, 1 priority 2
  void schedule_priority_01()
  {
    Solver s;
    vector<cspvar> d0, d1, l, u, f;
    vector<cspvar> x = s.newCSPVarArray(3, 5, 10);
    d0.assign(x.begin(), x.begin()+2);
    d1.assign(x.begin()+1, x.begin()+3);

    vector<int> wakes;
    action_schedule a1, a2, a3, a4, a5;
    post_test(s, 1, wakes,  x, l, u, f, a1, 2);
    post_test(s, 2, wakes, d0, l, u, f, a2, 1);
    post_test(s, 3, wakes, d1, l, u, f, a3, 1);
    post_test(s, 4, wakes, d0, l, u, f, a4, 0);
    post_test(s, 5, wakes, d1, l, u, f, a5, 0);

    post_eq(s, x[0], x[2], 0);

    int exp[] = { 4, 5, 2, 3, 1, -1 };
    s.newDecisionLevel();
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule_priority_01);

  // 2 priority 0, 2 priority 1, 1 priority 2, with some actions
  // thrown in
  void schedule_priority_02()
  {
    Solver s;
    vector<cspvar> d0, d1, l, u, f;
    vector<cspvar> x = s.newCSPVarArray(3, 5, 10);
    d0.assign(x.begin(), x.begin()+2);
    d1.assign(x.begin()+1, x.begin()+3);

    vector<int> wakes;
    action_schedule a1, a2, a3, a4, a5;
    a2[3].push_back( domevent( x[0], domevent::NEQ, 7 ) );
    a1[8].push_back( domevent( x[0], domevent::NEQ, 8 ) );

    post_test(s, 1, wakes,  x, l, u, f, a1, 2);
    post_test(s, 2, wakes, d0, l, u, f, a2, 1);
    post_test(s, 3, wakes, d1, l, u, f, a3, 1);
    post_test(s, 4, wakes, d0, l, u, f, a4, 0);
    post_test(s, 5, wakes, d1, l, u, f, a5, 0);

    post_eq(s, x[0], x[2], 0);

    int exp[] = { 4, 5, 2, 4, 5, 3, 2, 1, 4, 5, 2, 3, 1, -1 };
    s.newDecisionLevel();
    x[0].remove(s, 6, NO_REASON);
    assert( !s.propagate() );
    assert( compare_events(wakes, exp) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(schedule_priority_02);

  // only high priority
  void schedule_priority_03()
  {
    for(int p = 1; p <= MAX_PRIORITY; ++p) {
      Solver s;
      vector<cspvar> d, l, u, f;
      d = s.newCSPVarArray(3, 5, 10);

      vector<int> wakes;
      action_schedule a;
      post_test(s, 1, wakes,  d, l, u, f, a, p);

      int exp[] = { 1, -1 };
      s.newDecisionLevel();
      d[0].remove(s, 6, NO_REASON);
      assert( !s.propagate() );
      assert( compare_events(wakes, exp) );
      s.cancelUntil(0);
    }
  }
  REGISTER_TEST(schedule_priority_03);
}

void schedule_test()
{
  cerr << "propagator scheduling tests\n";
  the_test_container().run();
}
