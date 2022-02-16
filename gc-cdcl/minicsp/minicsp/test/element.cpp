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

namespace {

  // prune in I
  void element01()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(0, 9);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X);
    s.propagate();

    s.newDecisionLevel();
    X[5].setmax(s, 4, NO_REASON);
    s.propagate();
    assert( !I.indomain(s, 5) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(element01);

  // prune in R because I is assigned
  void element02()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(0, 9);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X);
    s.propagate();

    for(int i = 0; i != 10; ++i)
      X[i].remove(s, i+1, NO_REASON);
    s.propagate();
    for(int i = 4; i != 10; ++i) {
      s.newDecisionLevel();
      I.assign(s, i, NO_REASON);
      s.propagate();
      assert( !R.indomain(s, i+1) );
      s.cancelUntil(0);
    }
  }
  REGISTER_TEST(element02);

  // prune in R because the value has been removed from all X[i]
  // s.t. i in D(I)
  void element03()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(0, 9);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X);
    s.propagate();

    s.newDecisionLevel();
    for(int i = 0; i != 10; ++i)
      X[i].remove(s, 8, NO_REASON);
    s.propagate();
    assert(!R.indomain(s, 8));
    s.cancelUntil(0);
  }
  REGISTER_TEST(element03);

  // Prune in X after assigning I
  void element04()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(0, 9);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X);
    s.propagate();

    s.newDecisionLevel();
    I.assign(s, 4, NO_REASON);
    R.remove(s, 8, NO_REASON);
    s.propagate();
    assert( !X[4].indomain(s, 8) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(element04);

  // unsat
  void element05()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(0, 9);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X);
    s.propagate();

    s.newDecisionLevel();
    for(int i = 0; i != 9; ++i)
      X[i].remove(s, 8, NO_REASON);
    I.remove(s, 9, NO_REASON);
    R.assign(s, 8, NO_REASON);
    assert(s.propagate());
    s.cancelUntil(0);
  }
  REGISTER_TEST(element05);

  // I is a unary var
  void element06()
  {
    Solver s;
    cspvar R = s.newCSPVar(3, 7);
    cspvar I = s.newCSPVar(2, 2);
    vector<cspvar> X = s.newCSPVarArray(4, 1, 8);
    for(int i = 0; i != 4; ++i) {
      X[i].setmin(s, 1+2*i, NO_REASON);
      X[i].setmax(s, 2+2*i, NO_REASON);
    }
    post_element(s, R, I, X);

    assert( !s.propagate() );
    assert( R.min(s) == 5 );
    assert( R.max(s) == 6 );
  }
  REGISTER_TEST(element06);

  // Prune in X after assigning I, offset=1
  void element07()
  {
    Solver s;
    cspvar R = s.newCSPVar(5, 10);
    cspvar I = s.newCSPVar(1, 10);
    vector<cspvar> X = s.newCSPVarArray(10, 1, 10);
    post_element(s, R, I, X, 1);
    s.propagate();

    s.newDecisionLevel();
    I.assign(s, 5, NO_REASON);
    R.remove(s, 8, NO_REASON);
    s.propagate();
    assert( !X[4].indomain(s, 8) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(element07);

  // a bug that was uncovered by a flatzinc model
  void element08()
  {
    Solver s;
    cspvar R = s.newCSPVar(1, 1);
    cspvar I = s.newCSPVar(1, 2);
    vector<cspvar> X = s.newCSPVarArray(2, 0, 1);
    X[0].setmin(s, 1, NO_REASON);
    X[1].setmax(s, 0, NO_REASON);
    post_element(s, R, I, X, 1);
    s.propagate();
    assert( I.max(s) == 1 );
  }
  REGISTER_TEST(element08);

  // another bug from flatzinc: post correct clauses when one of the
  // values inside the bounds is not in X[]
  void element09()
  {
    Solver s;
    cspvar R = s.newCSPVar(1, 4);
    cspvar I = s.newCSPVar(1, 3);
    vector<cspvar> X;
    X.push_back(s.newCSPVar(1, 2));
    X.push_back(s.newCSPVar(1, 3));
    X.push_back(s.newCSPVar(1, 3));
    post_element(s, R, I, X, 1);
    assert( !s.propagate() );
  }
  REGISTER_TEST(element09);
}

void element_test()
{
  cerr << "element tests\n";

  the_test_container().run();
}
