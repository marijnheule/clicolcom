/*************************************************************************
minicsp

Copyright 2014 George Katsirelos

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
#include "test.hpp"

using namespace std;

namespace {
    void lex_leq01()
    {
        Solver s;
        vector<cspvar> x = s.newCSPVarArray(10, 0, 2);
        vector<cspvar> y = s.newCSPVarArray(10, 0, 2);
        post_lex_leq(s, x, y);
        assert( !s.propagate() );
        for(size_t i = 0; i != x.size(); ++i) {
            assert(x[i].min(s) == 0);
            assert(x[i].max(s) == 2);
            assert(y[i].min(s) == 0);
            assert(y[i].max(s) == 2);
        }

        s.newDecisionLevel();
        x[0].setmin(s, 1, NO_REASON);
        assert( !s.propagate() );
        assert( y[0].min(s) == 1 );
        s.newDecisionLevel();
        y[0].setmax(s, 1, NO_REASON);
        y[1].setmax(s, 1, NO_REASON);
        assert( !s.propagate() );
        assert( x[0].max(s) == 1 );
        assert( x[1].max(s) == 1 );
        s.newDecisionLevel();
        y[1].setmax(s, 0, NO_REASON);
        x[3].setmin(s, 1, NO_REASON);
        y[3].setmax(s, 0, NO_REASON);
        assert( !s.propagate() );
        assert( x[2].max(s) == 1 );
        assert( y[2].min(s) == 1 );
        s.cancelUntil(0);
    }
    REGISTER_TEST(lex_leq01);

    void lex_less01()
    {
        Solver s;
        vector<cspvar> x = s.newCSPVarArray(3, 0, 2);
        vector<cspvar> y = s.newCSPVarArray(3, 0, 2);
        post_lex_less(s, x, y);
        assert( !s.propagate() );
        for(size_t i = 0; i != x.size(); ++i) {
            assert(x[i].min(s) == 0);
            assert(x[i].max(s) == 2);
            assert(y[i].min(s) == 0);
            assert(y[i].max(s) == 2);
        }

        s.newDecisionLevel();
        x[0].setmin(s, 1, NO_REASON);
        y[0].setmax(s, 1, NO_REASON);
        y[1].setmax(s, 0, NO_REASON);
        assert( !s.propagate() );
        assert( y[0].min(s) == 1 );
        assert( x[0].max(s) == 1 );
        assert( y[2].min(s) == 1 );
        assert( x[2].max(s) == 1 );
        s.cancelUntil(0);
    }
    REGISTER_TEST(lex_less01);
}


void lex_test()
{
  cerr << "lexicographic ordering tests\n";

  the_test_container().run();
}
