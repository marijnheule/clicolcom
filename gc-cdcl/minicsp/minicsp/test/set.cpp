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
#include "minicsp/core/setcons.hpp"

namespace {
  void encoding01()
  {
    Solver s;
    s.newSetVar(1, 3);
    assert_num_solutions(s, 8);
  }
  REGISTER_TEST(encoding01);

  void setdiff01()
  {
    Solver s;
    setvar A = s.newSetVar(1, 7);
    setvar B = s.newSetVar(2, 6);
    setvar C = s.newSetVar(0, 5);
    post_setdiff(s, A, B, C);
    assert( !C.includes(s, 0) );
    assert( !A.includes(s, 7) );
  }
  REGISTER_TEST(setdiff01);

  // a couple of tests from flatzinc on which it initially barfed
  void setdiff02()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(1, 5);
    setvar C = s.newSetVar(2, 2);
    for(int i = 1; i <= 3; ++i)
      A.include(s, i, NO_REASON);
    B.include(s, 1, NO_REASON);
    B.include(s, 3, NO_REASON);
    B.include(s, 5, NO_REASON);
    B.exclude(s, 2, NO_REASON);
    B.exclude(s, 4, NO_REASON);
    post_setdiff(s, A, B, C);
  }
  REGISTER_TEST(setdiff02);

  void setdiff03()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(0, 1);
    setvar C = s.newSetVar(0, 4);
    post_setdiff(s, A, B, C);

    s.newDecisionLevel();
    A.include(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( C.includes(s, 3) );
    s.cancelUntil(0);

    assert_num_solutions(s, 32);
  }
  REGISTER_TEST(setdiff03);

  // non-overlapping universes: A and B
  void setdiff04()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(6, 9);
    setvar C = s.newSetVar(1, 9);
    post_setdiff(s, A, B, C);

    s.newDecisionLevel();
    A.include(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( C.includes(s, 1) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(setdiff04);

  // non-overlapping universes: everything
  void setdiff05()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(6, 9);
    setvar C = s.newSetVar(11, 15);
    post_setdiff(s, A, B, C);

    for(int i = 1; i <= 15; ++i) {
      assert( A.excludes(s, i) );
      assert( C.excludes(s, i) );
    }
  }
  REGISTER_TEST(setdiff05);

  // non-overlapping universes: C with A/B, reduces to A \subseteq B
  void setdiff06()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(2, 5);
    setvar C = s.newSetVar(11, 15);
    post_setdiff(s, A, B, C);

    assert( A.excludes(s, 1) );

    for(int i = 11; i <= 15; ++i) {
      assert( C.excludes(s, i) );
    }

    s.newDecisionLevel();
    A.include(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( B.includes(s, 3) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(setdiff06);

  void setsymdiff01()
  {
    Solver s;
    setvar A = s.newSetVar(1, 1);
    setvar B = s.newSetVar(2, 3);
    setvar C = s.newSetVar(1, 3);

    post_setsymdiff(s, A, B, C);

    A.include(s, 1, NO_REASON);
    B.include(s, 2, NO_REASON);
    C.include(s, 3, NO_REASON);

    assert_num_solutions(s, 1);
  }
  REGISTER_TEST(setsymdiff01);

  void setsymdiff02()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(2, 4);
    setvar C = s.newSetVar(1, 4);
    post_setsymdiff(s, A, B, C);

    s.newDecisionLevel();
    A.exclude(s, 1, NO_REASON);
    assert(!s.propagate());
    assert( C.excludes(s, 1) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    B.exclude(s, 4, NO_REASON);
    assert(!s.propagate());
    assert( C.excludes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 2, NO_REASON);
    B.exclude(s, 2, NO_REASON);
    A.exclude(s, 3, NO_REASON);
    B.include(s, 3, NO_REASON);
    assert(!s.propagate());
    assert( C.includes(s, 2) );
    assert( C.includes(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 2, NO_REASON);
    B.include(s, 2, NO_REASON);
    A.exclude(s, 3, NO_REASON);
    B.exclude(s, 3, NO_REASON);
    assert(!s.propagate());
    assert( C.excludes(s, 2) );
    assert( C.excludes(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 2, NO_REASON);
    B.include(s, 3, NO_REASON);
    assert(!s.propagate());
    assert( !C.includes(s, 2) );
    assert( !C.includes(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 2, NO_REASON);
    B.include(s, 3, NO_REASON);
    C.include(s, 2, NO_REASON);
    C.include(s, 3, NO_REASON);
    assert(!s.propagate());
    assert(A.excludes(s, 3));
    assert(B.excludes(s, 2));
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.exclude(s, 2, NO_REASON);
    B.exclude(s, 3, NO_REASON);
    C.include(s, 2, NO_REASON);
    C.include(s, 3, NO_REASON);
    assert(!s.propagate());
    assert(A.includes(s, 3));
    assert(B.includes(s, 2));
    s.cancelUntil(0);

    s.newDecisionLevel();
    C.include(s, 1, NO_REASON);
    C.include(s, 4, NO_REASON);
    assert(!s.propagate());
    assert(A.includes(s, 1));
    assert(B.includes(s, 4));
    s.cancelUntil(0);
  }
  REGISTER_TEST(setsymdiff02);

  // non-overlapping universes C and A/B
  void setsymdiff03()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(2, 4);
    setvar C = s.newSetVar(7, 10);
    post_setsymdiff(s, A, B, C);

    assert( A.excludes(s, 1) );
    assert( B.excludes(s, 4) );

    s.newDecisionLevel();
    A.include(s, 2, NO_REASON);
    assert(!s.propagate());
    assert( B.includes(s, 2) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    B.include(s, 2, NO_REASON);
    assert(!s.propagate());
    assert( A.includes(s, 2) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.exclude(s, 2, NO_REASON);
    assert(!s.propagate());
    assert( B.excludes(s, 2) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(setsymdiff03);

  void setsymdiff04()
  {
    Solver s;
    setvar A = s.newSetVar(1, 3);
    setvar B = s.newSetVar(4, 6);
    setvar C = s.newSetVar(7, 10);
    post_setsymdiff(s, A, B, C);

    for(int i = 1; i <= 6; ++i) {
      assert( A.excludes(s, i) );
      assert( B.excludes(s, i) );
    }
  }
  REGISTER_TEST(setsymdiff04);

  void seteq01()
  {
    Solver s;
    setvar A = s.newSetVar(0, 3);
    setvar B = s.newSetVar(1, 4);
    post_seteq(s, A, B);
    assert_num_solutions(s, 8);
  }
  REGISTER_TEST(seteq01);

  void seteq02()
  {
    Solver s;
    setvar A = s.newSetVar(0, 3);
    setvar B = s.newSetVar(1, 4);
    post_seteq(s, A, B);

    A.card(s).setmax(s, 2, NO_REASON);

    s.newDecisionLevel();
    A.include(s, 1, NO_REASON);
    assert(!s.propagate());
    assert( B.includes(s, 1) );

    s.newDecisionLevel();
    A.exclude(s, 2, NO_REASON);
    assert(!s.propagate());
    assert( B.excludes(s, 2));
    s.cancelUntil(0);

    assert_num_solutions(s, 7);
  }
  REGISTER_TEST(seteq02);

  void setneq01()
  {
    Solver s;
    setvar A = s.newSetVar(0,0);
    setvar B = s.newSetVar(0,0);
    post_setneq(s, A, B);

    s.newDecisionLevel();
    A.exclude(s, 0, NO_REASON);
    B.exclude(s, 0, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 0, NO_REASON);
    B.include(s, 0, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    assert_num_solutions(s, 2);
  }
  REGISTER_TEST(setneq01);

  void setneq02()
  {
    Solver s;
    setvar A = s.newSetVar(3,4);
    setvar B = s.newSetVar(2,3);
    post_setneq(s, A, B);

    s.newDecisionLevel();
    A.include(s, 3, NO_REASON);
    A.exclude(s, 4, NO_REASON);
    B.exclude(s, 2, NO_REASON);
    assert( !s.propagate() );
    assert( B.excludes(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.include(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( !B.excludes(s, 3) );

    s.newDecisionLevel();
    B.exclude(s, 2, NO_REASON);
    A.exclude(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( B.excludes(s, 3) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(setneq02);

  void setneq03()
  {
    Solver s;
    setvar A = s.newSetVar(3,5);
    setvar B = s.newSetVar(2,4);
    post_setneq(s, A, B);

    assert_num_solutions(s, 60);
  }
  REGISTER_TEST(setneq03);

  void setneq04()
  {
    Solver s;
    setvar a = s.newSetVar(1, 2);
    setvar b = s.newSetVar(4, 6);
    post_setneq(s, a, b);

    setvar c = s.newSetVar(9, 10);
    post_setneq(s, b, c);

    // just the fact that we are able to post these is enough
  }
  REGISTER_TEST(setneq04);

  void seteq_re01()
  {
    Solver s;
    setvar A = s.newSetVar(1, 2);
    setvar B = s.newSetVar(2, 3);
    Var r = s.newVar();
    post_seteq_re(s, A, B, Lit(r));

    s.newDecisionLevel();
    A.include(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( s.value(r) == l_False );
    s.cancelUntil(0);

    s.newDecisionLevel();
    B.include(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( s.value(r) == l_False );
    s.cancelUntil(0);

    s.newDecisionLevel();
    A.exclude(s, 1, NO_REASON);
    A.include(s, 2, NO_REASON);
    B.exclude(s, 3, NO_REASON);
    assert( !s.propagate() );

    s.newDecisionLevel();
    s.uncheckedEnqueue( Lit(r) );
    assert( !s.propagate() );
    assert( B.includes(s, 2) );
    s.cancelUntil(1);

    s.newDecisionLevel();
    s.uncheckedEnqueue( ~Lit(r) );
    assert( !s.propagate() );
    assert( B.excludes(s, 2) );
    s.cancelUntil(1);

    s.newDecisionLevel();
    B.include(s, 2, NO_REASON);
    assert( !s.propagate() );
    assert( s.value(r) == l_True );
    s.cancelUntil(1);

    s.newDecisionLevel();
    B.exclude(s, 2, NO_REASON);
    assert( !s.propagate() );
    assert( s.value(r) == l_False );
    s.cancelUntil(0);
  }
  REGISTER_TEST(seteq_re01);

  void seteq_re02()
  {
    Solver s;
    setvar a = s.newSetVar(1, 2);
    setvar b = s.newSetVar(5, 8);
    cspvar r = s.newCSPVar(0, 1);
    post_seteq_re(s, a, b, r);

    // the following cases have been tested already in 01, the
    // important thing is that we post the constraint in the first
    // place
    s.newDecisionLevel();
    a.include(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( r.max(s) == 0 );
    s.cancelUntil(0);

    s.newDecisionLevel();
    b.include(s, 5, NO_REASON);
    assert( !s.propagate() );
    assert( r.max(s) == 0 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(seteq_re02);

  void setin01()
  {
    Solver s;
    cspvar x = s.newCSPVar(1, 5);
    setvar a = s.newSetVar(2, 4);
    post_setin(s, x, a);
    assert( x.min(s) == 2 );
    assert( x.max(s) == 4 );

    s.newDecisionLevel();
    a.exclude(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( !x.indomain(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    x.remove(s, 4, NO_REASON);
    a.exclude(s, 2, NO_REASON);
    a.exclude(s, 3, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // 4 solutions for each value of x
    assert_num_solutions(s, 12);
  }
  REGISTER_TEST(setin01);

  void setin_re01()
  {
    Solver s;
    cspvar x = s.newCSPVar(1, 5);
    setvar a = s.newSetVar(2, 4);
    Var b = s.newVar();

    post_setin_re(s, x, a, Lit(b));

    // set b false
    s.newDecisionLevel();
    s.enqueue( ~Lit(b) );
    x.assign(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( a.excludes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    s.enqueue( ~Lit(b) );
    a.include(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( !x.indomain(s, 3) );
    s.cancelUntil(0);

    // set b true
    s.newDecisionLevel();
    s.enqueue( Lit(b) );
    a.exclude(s, 3, NO_REASON);
    assert( !s.propagate() );
    assert( !x.indomain(s, 3) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    s.enqueue( Lit(b) );
    x.remove(s, 4, NO_REASON);
    a.exclude(s, 2, NO_REASON);
    a.exclude(s, 3, NO_REASON);
    assert( s.propagate() );
    s.cancelUntil(0);

    // propagate to b
    s.newDecisionLevel();
    x.assign(s, 4, NO_REASON);
    a.exclude(s, 4, NO_REASON);
    assert(!s.propagate() );
    assert( s.value(b) == l_False );
    s.cancelUntil(0);

    s.newDecisionLevel();
    x.assign(s, 3, NO_REASON);
    a.include(s, 3, NO_REASON);
    assert(!s.propagate() );
    assert( s.value(b) == l_True );
    s.cancelUntil(0);
  }
  REGISTER_TEST(setin_re01);

  void set_intersect01()
  {
    Solver s;
    setvar a = s.newSetVar(1, 5);
    setvar b = s.newSetVar(3, 7);
    setvar c = s.newSetVar(1, 7);
    post_setintersect(s, a, b, c);

    assert(c.excludes(s, 1));
    assert(c.excludes(s, 2));
    assert(c.excludes(s, 6));
    assert(c.excludes(s, 7));

    s.newDecisionLevel();
    a.exclude(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( c.excludes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    c.include(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( a.includes(s, 4) );
    assert( b.includes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    c.exclude(s, 4, NO_REASON);
    a.include(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( b.excludes(s, 4) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_intersect01);

  void set_union01()
  {
    Solver s;
    setvar a = s.newSetVar(1, 5);
    setvar b = s.newSetVar(3, 7);
    setvar c = s.newSetVar(1, 7);
    post_setunion(s, a, b, c);

    s.newDecisionLevel();
    a.exclude(s, 4, NO_REASON);
    b.exclude(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( c.excludes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    c.include(s, 4, NO_REASON);
    b.exclude(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( a.includes(s, 4) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    c.exclude(s, 4, NO_REASON);
    assert( !s.propagate() );
    assert( a.excludes(s, 4) );
    assert( b.excludes(s, 4) );
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_union01);

  void set_subset01()
  {
    Solver s;
    setvar a = s.newSetVar(5, 12);
    setvar b = s.newSetVar(7, 11);
    post_setsubset(s, a, b);
    assert( a.excludes(s, 5) );
    assert( a.excludes(s, 6) );
    assert( a.excludes(s, 12) );

    s.newDecisionLevel();
    b.exclude(s, 8, NO_REASON);
    assert(!s.propagate());
    assert( a.excludes(s, 8) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 8, NO_REASON);
    assert(!s.propagate());
    assert( b.includes(s, 8) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 9, NO_REASON);
    b.exclude(s, 9, NO_REASON);
    assert(s.propagate());
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_subset01);

  void set_subseteq_re01()
  {
    Solver s;
    setvar a = s.newSetVar(5, 12);
    setvar b = s.newSetVar(7, 11);
    cspvar r = s.newCSPVar(0, 1);
    post_setsubseteq_re(s, a, b, r);

    s.newDecisionLevel();
    a.include(s, 6, NO_REASON);
    assert(!s.propagate());
    assert( r.max(s) == 0 );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 7, NO_REASON);
    b.exclude(s, 7, NO_REASON);
    assert(!s.propagate());
    assert( r.max(s) == 0 );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 7, NO_REASON);
    r.setmin(s, 1, NO_REASON);
    assert(!s.propagate());
    assert( b.includes(s, 7) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    for(int i = 5; i <= 12; ++i)
      if( i != 7 )
        a.exclude(s, i, NO_REASON);
    for(int i = 7; i <= 11; ++i)
      if( i != 7 )
        b.exclude(s, i, NO_REASON);
    a.include(s, 7, NO_REASON);
    r.setmax(s, 0, NO_REASON);
    assert(!s.propagate());
    assert( b.excludes(s, 7) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    for(int i = 5; i <= 12; ++i)
      if( i != 7 )
        a.exclude(s, i, NO_REASON);
    for(int i = 7; i <= 11; ++i)
      if( i != 7 )
        b.exclude(s, i, NO_REASON);
    b.include(s, 7, NO_REASON);
    r.setmax(s, 0, NO_REASON);
    assert(s.propagate());
    s.cancelUntil(0);

    s.newDecisionLevel();
    for(int i = 5; i <= 12; ++i)
      if( i != 7 )
        a.exclude(s, i, NO_REASON);
    for(int i = 7; i <= 11; ++i)
      if( i != 7 )
        b.exclude(s, i, NO_REASON);
    b.exclude(s, 7, NO_REASON);
    r.setmax(s, 0, NO_REASON);
    assert(!s.propagate());
    assert( a.includes(s, 7) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    for(int i = 5; i <= 12; ++i)
      if( i != 7 )
        a.exclude(s, i, NO_REASON);
    b.include(s, 7, NO_REASON);
    assert(!s.propagate());
    assert( r.min(s) == 1 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_subseteq_re01);

  void set_subseteq_re02()
  {
    Solver s;
    setvar a = s.newSetVar(1, 2);
    setvar b = s.newSetVar(4, 6);
    cspvar r = s.newCSPVar(0, 1);
    post_setsubseteq_re(s, a, b, r);

    s.newDecisionLevel();
    r.setmin(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( a.excludes(s, 1) );
    assert( a.excludes(s, 2) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( r.max(s) == 0 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_subseteq_re02);

  void set_subseteq_re03()
  {
    Solver s;
    setvar a = s.newSetVar(9, 10);
    setvar b = s.newSetVar(4, 6);
    cspvar r = s.newCSPVar(0, 1);
    post_setsubseteq_re(s, a, b, r);

    s.newDecisionLevel();
    r.setmin(s, 1, NO_REASON);
    assert( !s.propagate() );
    assert( a.excludes(s, 9) );
    assert( a.excludes(s, 10) );
    s.cancelUntil(0);

    s.newDecisionLevel();
    a.include(s, 9, NO_REASON);
    assert( !s.propagate() );
    assert( r.max(s) == 0 );
    s.cancelUntil(0);
  }
  REGISTER_TEST(set_subseteq_re03);
}


void set_test()
{
  std::cerr << "set tests\n";

  the_test_container().run();
}
