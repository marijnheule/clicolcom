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

#include <vector>
#include <iostream>

#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"
#include "test.hpp"

using namespace std;
using namespace regular;

namespace {
  void buildd(int da[][3], vector<transition>& d)
  {
    for(size_t i = 0; da[i][0] >= 0; ++i)
      d.push_back( transition(da[i][0], da[i][1], da[i][2]) );
  }

  // non-minimal DFA
  void regular01()
  {
    Solver s;
    int da[][3] = {
      {1, 0, 1},
      {1, 1, 2},
      {1, 2, 3},
      {2, 0, 2},
      {2, 1, 2},
      {2, 2, 2},
      {3, 0, 3},
      {3, 1, 3},
      {3, 2, 3},
      {-1, -1, -1}
    };
    vector<transition> d;
    set<int> f;
    buildd(da, d);
    f.insert(2);
    automaton a(d, 1, f);

    vector<cspvar> X = s.newCSPVarArray(5, 0, 3);
    post_regular(s, X, a, false);

    // symbol not in the language, fail
    s.newDecisionLevel();
    X[3].assign(s, 3, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // bad assignments, fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 2, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // non-failing assignments, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[2].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);

    // non-accepting complete assignment, fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 0, NO_REASON);
    X[3].assign(s, 0, NO_REASON);
    X[4].assign(s, 0, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // accepting partial assignment, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 1, NO_REASON);
    X[3].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);

    // accepting complete assignment, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 1, NO_REASON);
    X[3].assign(s, 2, NO_REASON);
    X[4].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);
  }
  REGISTER_TEST(regular01);

  // this tests for a bug in unfolding which occured when some early
  // state (i.e. state k with k < max_state_id) have no outgoing
  // transitions so do not appear in the transition table.
  //
  // The particular is example is exactly as in the test above with
  // states 2 and 3 renamed to each other and all the outgoing
  // transitions of the (absorbing, rejecting) state 2 removed. The
  // results should be exactly the same as above in terms of pruning
  // the Xs.
  //
  // This was uncovered while testing the xcsp frontend
  void regular02()
  {
    Solver s;
    int da[][3] = {
      {1, 0, 1},
      {1, 1, 3},
      {1, 2, 2},
      {3, 0, 3},
      {3, 1, 3},
      {3, 2, 3},
      {-1, -1, -1}
    };
    vector<transition> d;
    set<int> f;
    buildd(da, d);
    f.insert(3);
    automaton a(d, 1, f);

    vector<cspvar> X = s.newCSPVarArray(5, 0, 3);
    post_regular(s, X, a, false);

    // symbol not in the language, fail
    s.newDecisionLevel();
    X[3].assign(s, 3, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // bad assignments, fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 2, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // non-failing assignments, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[2].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);

    // non-accepting complete assignment, fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 0, NO_REASON);
    X[3].assign(s, 0, NO_REASON);
    X[4].assign(s, 0, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // accepting partial assignment, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 1, NO_REASON);
    X[3].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);

    // accepting complete assignment, don't fail
    s.newDecisionLevel();
    X[0].assign(s, 0, NO_REASON);
    X[1].assign(s, 0, NO_REASON);
    X[2].assign(s, 1, NO_REASON);
    X[3].assign(s, 2, NO_REASON);
    X[4].assign(s, 2, NO_REASON);
    assert( !s.propagate() );
    s.cancelUntil(0);
  }
  REGISTER_TEST(regular02);

  // some values unsupported at the root
  void regular03()
  {
      Solver s;
      int da[][3] = {
          {1, 1, 2},
          {1, 2, 5},
          {1, 0, 1},
          {2, 1, 2},
          {2, 2, 3},
          {2, 0, 2},
          {3, 1, 4},
          {3, 2, 3},
          {3, 0, 3},
          {4, 1, 4},
          {4, 2, 0},
          {4, 0, 4},
          {5, 1, 6},
          {5, 2, 5},
          {5, 0, 5},
          {6, 1, 6},
          {6, 2, 7},
          {6, 0, 6},
          {7, 1, 0},
          {7, 2, 7},
          {7, 0, 7},
          {-1, -1, -1}
      };

      vector<transition> d;
      set<int> f;
      buildd(da, d);
      for(int i = 1; i <= 7; ++i)
          f.insert(i);
      automaton a(d, 1, f);

      vector<cspvar> X = s.newCSPVarArray(8, 0, 2);
      X[0].assign(s, 1, NO_REASON);
      X[4].assign(s, 1, NO_REASON);
      X[6].assign(s, 2, NO_REASON);
      post_regular(s, X, a, true);
      assert(!s.propagate());

      for(int i = 1; i <= 3; ++i)
          assert( !X[i].indomain(s, 2) );
      assert(X[5].indomain(s, 2));
      assert(X[6].indomain(s, 2));
  }
  REGISTER_TEST(regular03);
}

void regular_test()
{
  cerr << "regular tests\n";
  the_test_container().run();
}
