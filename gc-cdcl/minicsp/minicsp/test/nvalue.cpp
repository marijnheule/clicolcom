/*************************************************************************
minicsp

Copyright 2011 George Katsirelos

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
  void nvalue01()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(3, 1, 5);
    cspvar N = s.newCSPVar(0, 3);

    post_atmostnvalue(s, x, N);
    assert( N.min(s) == 1 );

    s.newDecisionLevel();
    x[0].setmin(s, 3, NO_REASON);
    x[1].setmax(s, 2, NO_REASON);
    s.propagate();
    assert( N.min(s) == 2 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(nvalue01);

  // from nina's thesis
  void nvalue02()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(6, 1, 6);
    cspvar N = s.newCSPVar(1, 6);

    post_atmostnvalue(s, x, N);

    x[1].assign(s, 2, NO_REASON);
    x[2].assign(s, 2, NO_REASON);
    x[3].setmin(s, 2, NO_REASON);
    x[3].setmax(s, 4, NO_REASON);
    x[4].setmin(s, 4, NO_REASON);
    x[4].setmax(s, 5, NO_REASON);
    x[5].setmin(s, 4, NO_REASON);
    x[5].setmax(s, 5, NO_REASON);

    s.propagate();

    assert(N.min(s) == 2);

    s.newDecisionLevel();
    N.assign(s, 2, NO_REASON);
    s.propagate();
    assert(x[0].min(s) == 2);
    assert(x[0].max(s) == 5);
    s.cancelUntil(0);
  }
  REGISTER_TEST(nvalue02);

  // pruning min/max from the first/last/middle clique, pruning a var
  // than spans multiple cliques

  // before:
  //    1 2 3 4 5 6 7 8 9 A
  // X1   * *
  // X2         * *
  // X3               * *
  // X4       * * * *
  // X5   * * *
  // X6 * * *
  // X7             * * *
  // X8               * * *
  // X9       * * * * * * *
  // N = 3

  // after:
  //    1 2 3 4 5 6 7 8 9 A
  // X1   * *
  // X2         * *
  // X3               * *
  // X4         * *
  // X5   * *
  // X6   * *
  // X7               * *
  // X8               * *
  // X9         * * * * *
  // N = 3
  void nvalue03()
  {
    Solver s;
    s.debugclauses = 1;
    cspvar N = s.newCSPVar(1, 3);
    vector<cspvar> x = s.newCSPVarArray(9, 1, 10);

    post_atmostnvalue(s, x, N);

    s.newDecisionLevel();
    x[0].setmin(s, 2, NO_REASON);
    x[0].setmax(s, 3, NO_REASON);
    x[1].setmin(s, 5, NO_REASON);
    x[1].setmax(s, 6, NO_REASON);
    x[2].setmin(s, 8, NO_REASON);
    x[2].setmax(s, 9, NO_REASON);
    x[3].setmin(s, 4, NO_REASON);
    x[3].setmax(s, 7, NO_REASON);
    x[4].setmin(s, 2, NO_REASON);
    x[4].setmax(s, 4, NO_REASON);
    x[5].setmin(s, 1, NO_REASON);
    x[5].setmax(s, 3, NO_REASON);
    x[6].setmin(s, 7, NO_REASON);
    x[6].setmax(s, 9, NO_REASON);
    x[7].setmin(s, 8, NO_REASON);
    x[7].setmax(s, 10, NO_REASON);
    x[8].setmin(s, 4, NO_REASON);
    x[8].setmax(s, 10, NO_REASON);

    s.propagate();

    assert( x[3].min(s) == 5 );
    assert( x[3].max(s) == 6 );
    assert( x[4].min(s) == 2 );
    assert( x[4].max(s) == 3 );
    assert( x[5].min(s) == 2 );
    assert( x[5].max(s) == 3 );
    assert( x[6].min(s) == 8 );
    assert( x[6].max(s) == 9 );
    assert( x[7].min(s) == 8 );
    assert( x[7].max(s) == 9 );
    assert( x[8].min(s) == 5 );
    assert( x[8].max(s) == 9 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(nvalue03);
}

void nvalue_test()
{
  cerr << "nvalue tests\n";
  the_test_container().run();
}
