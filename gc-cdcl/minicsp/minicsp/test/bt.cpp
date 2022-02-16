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

#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"

using namespace std;
using namespace minicsp;

namespace {
  void test01()
  {
    Solver s;
    btptr p = s.alloc_backtrackable(sizeof(int));
    btptr p1 = s.alloc_backtrackable(sizeof(char));

    int & i = s.deref<int>(p);
    char & c = s.deref<char>(p1);

    i = 0;
    c = 1;

    s.newDecisionLevel();
    i = 1;
    c = 2;

    s.newDecisionLevel();
    i = 2;
    c = 3;

    s.cancelUntil(1);
    assert( i == 1 );
    assert( c == 2 );

    s.newDecisionLevel();
    i = 3;
    c = 4;

    s.newDecisionLevel();
    i = 4;
    c = 5;

    s.cancelUntil(2);
    assert(i == 3);
    assert(c == 4);

    s.cancelUntil(1);
    assert( i == 1 );
    assert( c == 2 );

    s.cancelUntil(0);
    assert( i == 0 );
    assert( c == 1 );
  }
}

void bt_test()
{
  cerr << "bt tests\n";

  cerr << "test 01 ... " << flush;
  test01();
  cerr << "OK" << endl;
}
