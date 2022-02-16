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

#include <iostream>
#include "test.hpp"
#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"

using std::vector;

namespace {
  /* a or -b or -c
     ast: -a, b, c
     exp: fail, clause (a, -b, -c)
  */
  void test01()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ 1, -1, -1 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, -1);

    // assignment
    s.enqueue(~Lit(v[0]), 0L);
    s.enqueue(Lit(v[1]), 0L);
    s.enqueue(Lit(v[2]), 0L);

    // clause
    vec<Lit> expected;
    expected.push( Lit(v[0]) );
    expected.push( ~Lit(v[1]) );
    expected.push( ~Lit(v[2]) );

    Clause *c = s.propagate();
    assert(c);
    assert_clause_exact(s, c, expected);
  }
  REGISTER_TEST(test01);

  /* -2*a + b + c >= 0
     ast: a, -b
     exp: fail, clause (-a, b)
  */
  void test02()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ -2, 1, 1 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(Lit(v[0]), 0L);
    s.enqueue(~Lit(v[1]), 0L);

    // clause
    vec<Lit> expected;
    expected.push( ~Lit(v[0]) );
    expected.push( Lit(v[1]) );

    Clause *c = s.propagate();
    assert(c);
    assert_clause_exact(s, c, expected);
  }
  REGISTER_TEST(test02);

  /* -2*a + b + c >= 0
     ast: a
     exp: no fail
  */
  void test03()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ -2, 1, 1 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(Lit(v[0]), 0L);

    Clause *c = s.propagate();
    assert(!c);
  }
  REGISTER_TEST(test03);

  /* -100*a + 75*b + 75*c >= 0
     ast: a, b
     exp: no fail
  */
  void test04()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ -100, 75, 75 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(Lit(v[0]), 0L);
    s.enqueue(Lit(v[1]), 0L);

    Clause *c = s.propagate();
    assert(!c);
  }
  REGISTER_TEST(test04);

  /* -100*a + 75*b + 75*c >= 0
     ast: b, c
     exp: no fail
  */
  void test05()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ -100, 75, 75 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(Lit(v[1]), 0L);
    s.enqueue(Lit(v[2]), 0L);

    Clause *c = s.propagate();
    assert(!c);
  }
  REGISTER_TEST(test05);

  /* +100*a - 75*b - 75*c >= 0
     ast: b, c
     exp: fail, clause (-b, -c)
  */
  void test06()
  {
    Solver s;
    vector<Var> v(3);
    for(int i = 0; i != 3; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ 100, -75, -75 };
    vector<int> w(weights, weights+3);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(Lit(v[1]), 0L);
    s.enqueue(Lit(v[2]), 0L);

    // clause
    vec<Lit> expected;
    expected.push( ~Lit(v[1]) );
    expected.push( ~Lit(v[2]) );

    Clause *c = s.propagate();
    assert(c);
    assert_clause_exact(s, c, expected);
  }
  REGISTER_TEST(test06);

  /* +100*a - 75*b - 75*c + 10*d >= 0
     ast: -d, b, c
     exp: fail, clause (-b, -c)
     tests clauses minimality
  */
  void test07()
  {
    Solver s;
    vector<Var> v(4);
    for(int i = 0; i != 4; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ 100, -75, -75, 10 };
    vector<int> w(weights, weights+4);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(~Lit(v[3]), 0L);
    s.enqueue(Lit(v[1]), 0L);
    s.enqueue(Lit(v[2]), 0L);

    // clause
    vec<Lit> expected;
    expected.push( ~Lit(v[1]) );
    expected.push( ~Lit(v[2]) );

    Clause *c = s.propagate();
    assert(c);
    assert_clause_exact(s, c, expected);
  }
  REGISTER_TEST(test07);

  /* +10*a - 75*b - 75*c + 100*d >= 0
     ast: -a, b, c
     exp: fail, clause (-b, -c)
     tests clauses minimality and _svars sorting
  */
  void test08()
  {
    Solver s;
    vector<Var> v(4);
    for(int i = 0; i != 4; ++i)
      v[i] = s.newVar();

    // constraint
    int weights[]={ 10, -75, -75, 100 };
    vector<int> w(weights, weights+4);
    post_pb(s, v, w, 0);

    // assignment
    s.enqueue(~Lit(v[0]), 0L);
    s.enqueue(Lit(v[1]), 0L);
    s.enqueue(Lit(v[2]), 0L);

    // clause
    vec<Lit> expected;
    expected.push( ~Lit(v[1]) );
    expected.push( ~Lit(v[2]) );

    Clause *c = s.propagate();
    assert(c);
    assert_clause_exact(s, c, expected);
  }
  REGISTER_TEST(test08);

  /* 5*x1+5*x2 >= 0 ==> b = 1
     exp: b=1
   */
  void test09()
  {
    // test disabled for now because the pb constraint does not do
    // propagation
#if 0
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(2, 0, 1);
    cspvar b = s.newCSPVar(0, 1);
    vector<int> w(2);
    w[0] = 5;
    w[1] = 5;
    post_pb_right_imp_re(s, x, w, 0, b);
#endif
  }
  REGISTER_TEST(test09);

  /* pbvar, unit weights, c = 0 */
  void pbvar01()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 0, 1);
    vector<int> w(5);
    for(int i = 0; i != 5; ++i) w[i] = 1;
    cspvar rhs = s.newCSPVar(-5, 10);
    post_pb(s, x, w, 0, rhs);
    assert( rhs.min(s) == 0 );
    assert( rhs.max(s) == 5 );

    s.newDecisionLevel();
    x[0].setmin(s, 1, NO_REASON);
    x[1].setmax(s, 0, NO_REASON);
    s.propagate();
    assert(rhs.min(s) == 1);
    assert(rhs.max(s) == 4);

    s.newDecisionLevel();
    rhs.setmax(s, 1, NO_REASON);
    assert( !s.propagate() );
    for(int i = 2; i != 5; ++i)
      assert(x[i].max(s) == 0);
    s.cancelUntil(0);
  }
  REGISTER_TEST(pbvar01);

  void pbvar02()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 0, 1);
    vector<int> w(5);
    for(int i = 0; i != 5; ++i) w[i] = 1;
    cspvar rhs = s.newCSPVar(-5, 10);
    post_pb(s, x, w, 0, rhs);
    assert( rhs.min(s) == 0 );
    assert( rhs.max(s) == 5 );

    s.newDecisionLevel();
    rhs.setmin(s, 2, NO_REASON);
    x[0].setmin(s, 1, NO_REASON);
    x[1].setmax(s, 0, NO_REASON);
    x[2].setmax(s, 0, NO_REASON);
    x[3].setmax(s, 0, NO_REASON);
    assert( !s.propagate() );
    assert(x[4].min(s) == 1);
    assert(rhs.max(s) == 2);
    s.cancelUntil(0);
  }
  REGISTER_TEST(pbvar02);

  void pbvar03()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 0, 1);
    vector<int> w(5);
    for(int i = 0; i != 5; ++i) w[i] = 1;
    w[4] = -1;
    cspvar rhs = s.newCSPVar(-5, 10);
    post_pb(s, x, w, 0, rhs);
    assert( rhs.min(s) == -1 );
    assert( rhs.max(s) == 4 );

    s.newDecisionLevel();
    rhs.setmin(s, 1, NO_REASON);
    x[0].setmin(s, 1, NO_REASON);
    x[1].setmax(s, 0, NO_REASON);
    x[2].setmax(s, 0, NO_REASON);
    x[3].setmax(s, 0, NO_REASON);
    assert( !s.propagate() );
    assert(x[4].max(s) == 0);
    assert(rhs.max(s) == 1);
    s.cancelUntil(0);
  }
  REGISTER_TEST(pbvar03);

  // this is meant to test the clause produced, so is meaningful with
  // EXPENSIVE_INVARIANTS turned on.
  void pbvar04()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(3, 0, 1);
    vector<int> w(3);
    for(int i = 0; i != 3; ++i) w[i] = 1;
    cspvar rhs = s.newCSPVar(0, 3);
    post_pb(s, x, w, 0, rhs);

    s.newDecisionLevel();
    rhs.remove(s, 1, NO_REASON);
    x[1].setmax(s, 0, NO_REASON);
    x[2].setmax(s, 0, NO_REASON);
    assert( !s.propagate() );
    assert(x[0].max(s) == 0);
    assert(rhs.max(s) == 0);
    s.cancelUntil(0);
  }
  REGISTER_TEST(pbvar04);

  // test that pruning the rhs does not put any lit_undefs in the clause
  void pbvar05()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(3, 0, 1);
    vector<int> w(3);
    for(int i = 0; i != 3; ++i) w[i] = 3;
    w[1] = 0;
    cspvar rhs = s.newCSPVar(0, 3);
    post_pb(s, x, w, 0, rhs);

    s.newDecisionLevel();
    x[0].setmin(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert(rhs.min(s) == 3);
    assert(rhs.max(s) == 3);
    s.cancelUntil(0);
  }
  REGISTER_TEST(pbvar05);
}

void pb_test()
{
  std::cerr << "pb tests\n";

  the_test_container().run();
}
