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
  void test01()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, -5, 5);
    vector<int> w(5);
    w[0] = 5;
    w[1] = -4;
    w[2] = 3;
    w[3] = -2;
    w[4] = 1;
    post_lin_leq(s, x, w, 0);
    s.propagate();
  }
  REGISTER_TEST(test01);

  /* sum x1 .. x5 >= 24
     x1 .. x5 \in [3,5]

     result x1..x5 >= 4
   */
  void test02()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 3, 5);
    vector<int> w(5);
    w[0] = -1;
    w[1] = -1;
    w[2] = -1;
    w[3] = -1;
    w[4] = -1;
    post_lin_leq(s, x, w, 24);
    s.propagate();
     for(int i = 0; i != 5; ++i) {
      assert( x[i].min(s) == 4 );
    }
  }
  REGISTER_TEST(test02);

  /* sum x1 .. x5 <= 20
     x1 .. x5 \in [3, 9]

     result x1..x5<=8
   */
  void test03()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 3, 9);
    vector<int> w(5);
    w[0] = 1;
    w[1] = 1;
    w[2] = 1;
    w[3] = 1;
    w[4] = 1;
    post_lin_leq(s, x, w, -20);
    s.propagate();
    for(int i = 0; i != 5; ++i) {
      assert( x[i].max(s) == 8 );
    }
  }
  REGISTER_TEST(test03);

  /* sum x1...x4 - x5 <= 20
     x1 .. x5 \in [7,10]
     x6 \in [2, 10]

     result x5 >= 8
     x1..x4 <= 9
   */
  void test04()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 2, 10);
    vector<int> w(5);
    w[0] = 1;
    w[1] = 1;
    w[2] = 1;
    w[3] = 1;
    w[4] = -1;
    for(int i = 0; i != 4; ++i)
      x[i].setmin(s, 7, NO_REASON);
    post_lin_leq(s, x, w, -20);
    s.propagate();
    for(int i = 0; i != 4; ++i) {
      assert( x[i].max(s) == 9 );
    }
    assert( x[4].min(s) == 8 );
  }
  REGISTER_TEST(test04);

  /* sum x1 .. x4 <= 20 implies b
     x1 .. x4 \in [3, 5]

     result: b true
   */
  void test05()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 3, 5);
    vector<int> w(4);
    cspvar b = s.newCSPVar(0, 1);
    w[0] = 1;
    w[1] = 1;
    w[2] = 1;
    w[3] = 1;
    post_lin_leq_right_imp_re(s, x, w, -20, b);
    s.propagate();
    assert( b.min(s) == 1 );
  }
  REGISTER_TEST(test05);

  /* sum x1 .. x4 <= 20 implies b
     x1 .. x4 \in [1,6]

     branch 1:
     set all xi <= 5, result: b true

     branch 2:
     set b false, result: all xi >= 3
  */
  void test06()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 1, 6);
    vector<int> w(4);
    cspvar b = s.newCSPVar(0, 1);
    w[0] = 1;
    w[1] = 1;
    w[2] = 1;
    w[3] = 1;
    post_lin_leq_right_imp_re(s, x, w, -20, b);
    Clause *confl = s.propagate();
    assert(!confl);

    s.newDecisionLevel();
    for(int i = 0; i != 4; ++i)
      x[i].setmax(s, 5, NO_REASON);
    confl = s.propagate();
    assert(!confl);
    assert(b.min(s) == 1);

    s.cancelUntil(0);
    s.newDecisionLevel();
    b.setmax(s, 0, NO_REASON);
    confl = s.propagate();
    assert(!confl);
    for(int i = 0; i != 4; ++i)
      assert( x[i].min(s) == 3 );
  }
  REGISTER_TEST(test06);

  /* b --> sum x1 .. x4 <= 20
     x1 .. x4 \in 6..10

     result: b false
   */
  void test07()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 6, 10);
    vector<int> w(4);
    cspvar b = s.newCSPVar(0, 1);
    w[0] = 1;
    w[1] = 1;
    w[2] = 1;
    w[3] = 1;
    post_lin_leq_left_imp_re(s, x, w, -20, b);
    s.propagate();
    assert(b.max(s) == 0);
  }
  REGISTER_TEST(test07);

  /* b --> sum -x1 .. -x4 <= -20
     x1 .. x4 \in [1,6]

     branch 1:
     set all xi <= 4, result: b false

     branch 2:
     set b true, result: all xi >= 2
   */
  void test08()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 1, 6);
    vector<int> w(4);
    cspvar b = s.newCSPVar(0, 1);
    w[0] = -1;
    w[1] = -1;
    w[2] = -1;
    w[3] = -1;
    post_lin_leq_left_imp_re(s, x, w, 20, b);
    Clause *confl = s.propagate();
    assert(!confl);

    s.newDecisionLevel();
    for(int i = 0; i != 4; ++i)
      x[i].setmax(s, 4, NO_REASON);
    s.propagate();
    assert(b.max(s) == 0);

    s.cancelUntil(0);
    s.newDecisionLevel();
    b.assign(s, 1, NO_REASON);
    s.propagate();
    for(int i = 0; i != 4; ++i)
      assert( x[i].min(s) == 2 );
  }
  REGISTER_TEST(test08);

  /* not-so-black-box. this tests that the clause produced from
     lin_leq contains the right literals, even though x[1].r_geq(1)
     and x[2].r_geq(1) are both lit_Undef. This would mean that
     x[3].e_leq(2) does not go to _ps[3] but to _ps[1] in the original
     buggy implementation.
  */
  void test09()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 1, 20);
    vector<int> c(4);
    c[0] = 1000;
    c[1] = 1;
    c[2] = 1;
    c[3] = 1000;
    post_lin_leq(s, x, c, -5040);
    x[0].setmin(s, 3, NO_REASON);
    x[3].setmin(s, 2, NO_REASON);
    Clause *confl = s.propagate();
    assert(!confl);
    assert(x[0].max(s) == 3);
    assert(x[3].max(s) == 2);
  }
  REGISTER_TEST(test09);

  /* sum x1..x4 <= 18 ===> b
     x1..x4 in [1..5]
     b = 0 *while posting*

     result: x1..x4 >= 4
   */
  void test10()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 1, 5);
    cspvar b = s.newCSPVar(0,1);
    vector<int> c(4);
    c[0] = 1;
    c[1] = 1;
    c[2] = 1;
    c[3] = 1;
    b.setmax(s, 0, NO_REASON);
    post_lin_leq_right_imp_re(s, x, c, -18, b);
    assert(x[0].min(s) == 4);
  }
  REGISTER_TEST(test10);

  /* b ===> sum x1..x4 <= 18
     x1..x4 in [4, 7]
     b = 1 *while posting*

     result x1..x4 <= 6
   */
  void test11()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(4, 4, 7);
    cspvar b = s.newCSPVar(0,1);
    vector<int> c(4);
    c[0] = 1;
    c[1] = 1;
    c[2] = 1;
    c[3] = 1;
    b.setmin(s, 1, NO_REASON);
    post_lin_leq_left_imp_re(s, x, c, -18, b);
    assert(x[0].max(s) == 6);
  }
  REGISTER_TEST(test11);

  /* sum x1..x5 <= 25
     x1 .. x4 in [4, 10]
     x5 = 5

     result: x1..x4 <= 8

     this is supposed to test that the correct constraint is posted
     after simplification when a variable is fixed. it's not clear how
     to test this other than making sure the correct propagation
     happens.
   */
  void test12()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(5, 4, 10);
    x[4].assign(s, 5, NO_REASON);
    vector<int> c(5);
    for(int i = 0; i != 5; ++i) c[i] = 1;
    post_lin_leq(s, x, c, -25);
  }
  REGISTER_TEST(test12);

  // unary constraint
  void test13()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(1, 5, 10);
    vector<int> c(1); c[0] = 1;
    post_lin_eq(s, x, c, -7);
    s.propagate();
    assert(x[0].max(s) == 7);
  }
  REGISTER_TEST(test13);

  /* SEND+MORE=MONEY */
  void test_money()
  {
    Solver s;
    vector<cspvar> x = s.newCSPVarArray(8, 0, 9);
    cspvar S = x[0]; s.setCSPVarName(S, "S");
    cspvar E = x[1]; s.setCSPVarName(E, "E");
    cspvar N = x[2]; s.setCSPVarName(N, "N");
    cspvar D = x[3]; s.setCSPVarName(D, "D");
    cspvar M = x[4]; s.setCSPVarName(M, "M");
    cspvar O = x[5]; s.setCSPVarName(O, "O");
    cspvar R = x[6]; s.setCSPVarName(R, "R");
    cspvar Y = x[7]; s.setCSPVarName(Y, "Y");

    vector<cspvar> v;
    vector<int> c;
    v.push_back(S); c.push_back(1000);
    v.push_back(E); c.push_back(100);
    v.push_back(N); c.push_back(10);
    v.push_back(D); c.push_back(1);
    v.push_back(M); c.push_back(1000);
    v.push_back(O); c.push_back(100);
    v.push_back(R); c.push_back(10);
    v.push_back(E); c.push_back(1);
    v.push_back(M); c.push_back(-10000);
    v.push_back(O); c.push_back(-1000);
    v.push_back(N); c.push_back(-100);
    v.push_back(E); c.push_back(-10);
    v.push_back(Y); c.push_back(-1);

    vector<int> c1(c);
    for(size_t i = 0; i != c1.size(); ++i) c1[i] = -c1[i];

    post_lin_leq(s, v, c, 0);
    post_lin_leq(s, v, c1, 0);
    post_alldiff(s, x);

    s.solve();
    cout << s.conflicts << " conflicts.\n";
    for(int i = 0; i != s.cspmodel.size(); ++i)
      assert(s.cspmodel[i].first == s.cspmodel[i].second);
    cout << "SEND+MORE=MONEY\n";
    vec< pair<int, int> >& m = s.cspmodel;
    cout << m[0].first << m[1].first << m[2].first << m[3].first << "+"
         << m[4].first << m[5].first << m[6].first << m[1].first << "="
         << m[4].first
         << m[5].first << m[2].first << m[1].first << m[7].first
         << "\n";
    int send = 1000*m[0].first + 100*m[1].first + 10*m[2].first
      + m[3].first;
    int more = 1000*m[4].first + 100*m[5].first + 10*m[6].first
      + m[1].first;
    int money = 10000*m[4].first +
      1000*m[5].first + 100*m[2].first + 10*m[1].first + m[7].first;
    assert(send+more==money);
  }
  REGISTER_TEST(test_money);
}

void lin_test()
{
  cerr << "lin tests\n";

  the_test_container().run();
}
