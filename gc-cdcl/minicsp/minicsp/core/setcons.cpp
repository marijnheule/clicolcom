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

#include "solver.hpp"
#include "setcons.hpp"
#include "cons.hpp"
#include <algorithm>

namespace minicsp {

void post_setdiff(Solver &s, setvar a, setvar b, setvar c)
{
  vec<Lit> ps;
  for(int i = c.umin(s); i < a.umin(s); ++i)
    c.exclude(s, i, NO_REASON);
  for(int i = a.umax(s)+1; i <= c.umax(s); ++i)
    c.exclude(s, i, NO_REASON);
  for(int i = a.umin(s); i <= a.umax(s); ++i) {
    if( i < b.umin(s) || i > b.umax(s) ) {
      // a.ini(i) <=> c.ini(i)
      if( i >= c.umin(s) && i <= c.umax(s) ) {
        ps.clear();
        ps.push( ~Lit( c.ini(s, i) ) );
        ps.push( Lit( a.ini(s, i) ) );
        s.addClause(ps);

        ps.clear();
        pushifdef( ps, Lit( c.ini(s, i) ) );
        ps.push( ~Lit( a.ini(s, i) ) );
        s.addClause(ps);
      } else {
        ps.clear();
        ps.push( ~Lit( a.ini(s, i) ) );
        s.addClause(ps);
      }
    } else {
      // a.ini(i) /\ !b.ini(i) => c.ini(i)
      ps.clear();
      pushifdef( ps, Lit( c.ini(s, i) ) );
      ps.push( Lit( b.ini(s, i) ) );
      ps.push( ~Lit( a.ini(s, i) ) );
      s.addClause(ps);

      // c.ini(i) => a.ini(i) /\ !b.ini(i))
      if( i < c.umin(s) || i > c.umax(s) )
        continue;

      ps.clear();
      ps.push( ~Lit( c.ini(s, i) ) );
      ps.push( Lit( a.ini(s, i) ) );
      s.addClause(ps);

      ps.clear();
      ps.push( ~Lit( c.ini(s, i) ) );
      ps.push( ~Lit( b.ini(s, i) ) );
      s.addClause(ps);
    }
  }
  post_leq(s, c.card(s), a.card(s), 0);
}

void post_setsymdiff(Solver& s, setvar a, setvar b, setvar c)
{
  vec<Lit> ps;
  for(int i = c.umin(s); i <= c.umax(s); ++i) {
    if( ( i < a.umin(s) || i > a.umax(s) ) &&
        ( i < b.umin(s) || i > b.umax(s) ) )
      c.exclude(s, i, NO_REASON);
  }
  for(int i = a.umin(s); i <= a.umax(s); ++i) {
    if( i < b.umin(s) || i > b.umax(s) ) {
      // a.ini(i) <=> c.ini(i)
      if( i >= c.umin(s) && i <= c.umax(s) ) {
        ps.clear();
        ps.push( ~Lit( c.ini(s, i) ) );
        ps.push( Lit( a.ini(s, i) ) );
        s.addClause(ps);

        ps.clear();
        ps.push( Lit( c.ini(s, i) ) );
        ps.push( ~Lit( a.ini(s, i) ) );
        s.addClause(ps);
      } else {
        ps.clear();
        ps.push( ~Lit( a.ini(s, i) ) );
        s.addClause(ps);
      }
    } else {
      // a.ini(i) /\ !b.ini(i) => c.ini(i)
      ps.clear();
      pushifdef( ps, Lit( c.ini(s, i) ) );
      ps.push( Lit( b.ini(s, i) ) );
      ps.push( ~Lit( a.ini(s, i) ) );
      s.addClause(ps);

      // !a.ini(i) /\ b.ini(i) => c.ini(i)
      ps.clear();
      pushifdef( ps, Lit( c.ini(s, i) ) );
      ps.push( ~Lit( b.ini(s, i) ) );
      ps.push( Lit( a.ini(s, i) ) );
      s.addClause(ps);

      if( i < c.umin(s) || i > c.umax(s) )
        continue;

      // c.ini(i) => a.ini(i) \/ b.ini(i))
      // c.ini(i) => (a.ini(i) => !b.ini(i))
      ps.clear();
      ps.push( ~Lit( c.ini(s, i) ) );
      ps.push( Lit( a.ini(s, i) ) );
      ps.push( Lit( b.ini(s, i) ) );
      s.addClause(ps);

      ps.clear();
      ps.push( ~Lit( c.ini(s, i) ) );
      ps.push( ~Lit( a.ini(s, i) ) );
      ps.push( ~Lit( b.ini(s, i) ) );
      s.addClause(ps);
    }
  }

  for(int i = b.umin(s); i <= b.umax(s); ++i) {
    if( i < a.umin(s) || i > a.umax(s) ) {
      // b.ini(i) <=> c.ini(i)
      if( i >= c.umin(s) && i <= c.umax(s) ) {
        ps.clear();
        ps.push( ~Lit( c.ini(s, i) ) );
        ps.push( Lit( b.ini(s, i) ) );
        s.addClause(ps);

        ps.clear();
        pushifdef( ps, Lit( c.ini(s, i) ) );
        ps.push( ~Lit( b.ini(s, i) ) );
        s.addClause(ps);
      } else {
        ps.clear();
        ps.push( ~Lit( b.ini(s, i) ) );
        s.addClause(ps);
      }
    } // else has been handled in the loop iterating over the universe
      // of a
  }
}

void post_seteq(Solver &s, setvar a, setvar b)
{
  post_eq(s, a.card(s), b.card(s), 0);

  for(int i = a.umin(s); i < b.umin(s); ++i)
    a.exclude(s, i, NO_REASON);
  for(int i = b.umax(s)+1; i <= a.umax(s); ++i)
    a.exclude(s, i, NO_REASON);
  for(int i = b.umin(s); i < a.umin(s); ++i)
    b.exclude(s, i, NO_REASON);
  for(int i = a.umax(s)+1; i <= b.umax(s); ++i)
    b.exclude(s, i, NO_REASON);

  vec<Lit> ps;
  for(int i = std::max(a.umin(s), b.umin(s));
      i <= std::min(a.umax(s), b.umax(s)); ++i) {
    ps.clear();
    ps.push( ~Lit( a.ini(s, i) ) );
    ps.push( Lit( b.ini(s, i) ) );
    s.addClause(ps);

    ps.clear();
    ps.push( Lit( a.ini(s, i) ) );
    ps.push( ~Lit( b.ini(s, i) ) );
    s.addClause(ps);
  }
}

namespace setneq {
  void post_unique_used(Solver &s, setvar a, setvar b,
                        vec<Lit>& ps1, vec<Lit>& ps2,
                        Var &onlya)
  {
    if( a.umin(s) < b.umin(s) || a.umax(s) > b.umax(s) ) {
      onlya = s.newVar();
      ps1.push( Lit(onlya) );
      ps1.growTo(2);
      ps2.push( ~Lit(onlya) );
    }

    for(int i = a.umin(s); i < std::min(a.umax(s)+1, b.umin(s)); ++i) {
      ps1.growTo(2);
      ps1[0] = Lit(onlya);
      ps1[1] =  ~Lit( a.ini(s, i) );
      s.addClause(ps1);
      ps2.push( Lit( a.ini(s, i) ) );
    }
    for(int i = std::max(a.umin(s), b.umax(s)+1); i <= a.umax(s); ++i) {
      ps1.growTo(2);
      ps1[0] = Lit(onlya);
      ps1[1] =  ~Lit( a.ini(s, i) );
      s.addClause(ps1);
      ps2.push( Lit( a.ini(s, i) ) );
    }
    if( ps2.size() > 0 )
      s.addClause(ps2);
  }
}

/* We introduce a variable y_i for each element i in both universes
 * such that y_i <=> a.ini(i) != b.ini(i) and variables onlya, onlyb
 * for elements that are only in a's (b's) universe. Finally, a big
 * clause (onlya, onlyb, y_1, ..., y_d) (assuming the shared universe
 * is 1..d)
 *
 * FIXME: if we find an instance that really wants a lot of setneq
 * constraints, check if it is better to post an automaton that
 * asserts at least one of onlya, onlyb, y_1, ..., y_d is
 * true. Propagation on it would be more local...
 */
void post_setneq(Solver &s, setvar a, setvar b)
{
  vec<Lit> ps1, ps2;

  Var onlya = var_Undef, onlyb = var_Undef;

  setneq::post_unique_used(s, a, b, ps1, ps2, onlya);

  ps1.clear();
  ps2.clear();

  setneq::post_unique_used(s, b, a, ps1, ps2, onlyb);

  ps1.clear();
  ps2.clear();

  if( onlya != var_Undef ) ps2.push( Lit(onlya) );
  if( onlyb != var_Undef ) ps2.push( Lit(onlyb) );

  for(int i = std::max(a.umin(s), b.umin(s)),
        iend = std::min(a.umax(s), b.umax(s))+1;
      i < iend; ++i) {
    Var y = s.newVar();
    ps2.push( Lit(y) );

    // y_i => (ai => ~bi)
    ps1.growTo(3);
    ps1[0] = ~Lit(y);
    ps1[1] = ~Lit( a.ini(s, i) );
    ps1[2] = ~Lit( b.ini(s, i) );
    s.addClause(ps1);

    // y_i => (~ai => bi)
    ps1.growTo(3);
    ps1[0] = ~Lit(y);
    ps1[1] = Lit( a.ini(s, i) );
    ps1[2] = Lit( b.ini(s, i) );
    s.addClause(ps1);

    // ~ai /\ bi => y_i
    ps1.growTo(3);
    ps1[0] = Lit(y);
    ps1[1] = Lit( a.ini(s, i) );
    ps1[2] = ~Lit( b.ini(s, i) );
    s.addClause(ps1);

    // ~bi /\ ai => y_i
    ps1.growTo(3);
    ps1[0] = Lit(y);
    ps1[1] = ~Lit( a.ini(s, i) );
    ps1[2] = Lit( b.ini(s, i) );
    s.addClause(ps1);
  }
  s.addClause(ps2);
}

/* Implemented the same way as setneq, only with the appropriate p =>
 * ... and ~p => ... added to the clauses.
 *
 * FIXME: make post_setneq forward to this if it is not too wasteful?
 */

void post_seteq_re(Solver &s, setvar a, setvar b, Lit p)
{
  vec<Lit> ps1, ps2;

  Var onlya = var_Undef, onlyb = var_Undef;

  setneq::post_unique_used(s, a, b, ps1, ps2, onlya);

  ps1.clear();
  ps2.clear();

  setneq::post_unique_used(s, b, a, ps1, ps2, onlyb);

  ps1.clear();
  ps2.clear();

  if( onlya != var_Undef ) {
    ps2.push( Lit(onlya) );
    ps1.growTo(2);
    ps1[0] = ~Lit(onlya);
    ps1[1] = ~p;
    s.addClause(ps1);
  }
  if( onlyb != var_Undef ) {
    ps2.push( Lit(onlyb) );
    ps1.growTo(2);
    ps1[0] = ~Lit(onlyb);
    ps1[1] = ~p;
    s.addClause(ps1);
  }

  for(int i = std::max(a.umin(s), b.umin(s)),
        iend = std::min(a.umax(s), b.umax(s))+1;
      i < iend; ++i) {
    Var y = s.newVar();
    ps2.push( Lit(y) );

    ps1.clear();
    ps1.growTo(2);
    ps1[0] = ~Lit(y);
    ps1[1] = ~p;
    s.addClause(ps1);

    // y_i => (ai => ~bi)
    ps1.growTo(3);
    ps1[0] = ~Lit(y);
    ps1[1] = ~Lit( a.ini(s, i) );
    ps1[2] = ~Lit( b.ini(s, i) );
    s.addClause(ps1);

    // y_i => (~ai => bi)
    ps1.growTo(3);
    ps1[0] = ~Lit(y);
    ps1[1] = Lit( a.ini(s, i) );
    ps1[2] = Lit( b.ini(s, i) );
    s.addClause(ps1);

    // ~ai /\ bi => y_i
    ps1.growTo(3);
    ps1[0] = Lit(y);
    ps1[1] = Lit( a.ini(s, i) );
    ps1[2] = ~Lit( b.ini(s, i) );
    s.addClause(ps1);

    // ~bi /\ ai => y_i
    ps1.growTo(3);
    ps1[0] = Lit(y);
    ps1[1] = ~Lit( a.ini(s, i) );
    ps1[2] = Lit( b.ini(s, i) );
    s.addClause(ps1);
  }
  ps2.push( p );
  s.addClause(ps2);
}

void post_seteq_re(Solver &s, setvar a, setvar b, cspvar r)
{
  assert(r.min(s) >= 0 && r.max(s) <= 1);
  if( r.min(s) == 1 ) {
    post_seteq(s, a, b);
    return;
  }
  if( r.max(s) == 0 ) {
    post_setneq(s, a, b);
    return;
  }
  post_seteq_re(s, a, b, Lit(r.eqi(s, 1)));
}


void post_setneq_re(Solver &s, setvar a, setvar b, Lit p)
{
  post_seteq_re(s, a, b, ~p);
}

void post_setneq_re(Solver &s, setvar a, setvar b, cspvar r)
{
  assert(r.min(s) >= 0 && r.max(s) <= 1);
  if( r.min(s) == 1 ) {
    post_setneq(s, a, b);
    return;
  }
  if( r.max(s) == 0 ) {
    post_seteq(s, a, b);
    return;
  }
  post_seteq_re(s, a, b, Lit(r.eqi(s, 0)));
}

void post_setin(Solver &s, cspvar x, setvar a)
{
  x.setmin(s, a.umin(s), NO_REASON);
  x.setmax(s, a.umax(s), NO_REASON);

  vec<Lit> ps;
  for(int i = x.min(s), iend = x.max(s)+1; i != iend; ++i) {
    ps.growTo(2);
    ps[0] = ~Lit( x.eqi(s, i) );
    ps[1] = Lit( a.ini(s, i) );
    s.addClause(ps);
  }
}

void post_setin_re(Solver &s, cspvar x, setvar a, Lit p)
{
  vec<Lit> ps;
  if( x.min(s) < a.umin(s) ) {
    ps.growTo(2);
    if( x.leqi(s, a.umin(s)-1) != var_Undef ) {
      ps[0] = ~Lit( x.leqi(s, a.umin(s)-1) );
      ps[1] = ~p;
      s.addClause(ps);
    } else {
      // x.max(s) < a.umin(s)
      s.enqueue(~p);
      return;
    }
  }

  if( x.max(s) > a.umax(s) ) {
    ps.growTo(2);
    if( x.leqi(s, a.umax(s)) != var_Undef ) {
      ps[0] = Lit( x.leqi(s, a.umax(s)) );
      ps[1] = ~p;
      s.addClause(ps);
    } else {
      // x.min(s) > a.umax(s)
      s.enqueue(~p);
      return;
    }
  }

  for(int i = std::max( x.min(s), a.umin(s) ),
        iend = std::min( x.max(s), a.umax(s) )+1;
      i != iend; ++i) {
    ps.growTo(3);
    ps[0] = ~p;
    ps[1] = ~Lit( x.eqi(s, i) );
    ps[2] = Lit( a.ini(s, i) );
    s.addClause(ps);

    ps.growTo(3);
    ps[0] = p;
    ps[1] = ~Lit( x.eqi(s, i) );
    ps[2] = ~Lit( a.ini(s, i) );
    s.addClause(ps);
  }
}

void post_setin_re(Solver &s, cspvar x, setvar a, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1);
  if( b.max(s) == 0 )
    post_setin_re(s, x, a, ~Lit( b.eqi(s, 0) ));
  else
    post_setin_re(s, x, a, Lit( b.eqi(s, 1) ));
}

void post_setintersect(Solver &s, setvar a, setvar b, setvar c)
{
  for(int i = c.umin(s); i < std::max(a.umin(s), b.umin(s)); ++i)
    c.exclude(s, i, NO_REASON);
  for(int i = std::min(a.umax(s), b.umax(s))+1; i <= c.umax(s); ++i)
    c.exclude(s, i, NO_REASON);

  vec<Lit> ps;
  for(int i = c.umin(s); i <= c.umax(s); ++i) {
    if( i < a.umin(s) || i < b.umin(s) ||
        i > a.umax(s) || i > b.umax(s) )
      continue;

    ps.clear();

    ps.growTo(2);
    ps[0] = Lit(a.ini(s, i));
    ps[1] = ~Lit(c.ini(s, i));
    s.addClause(ps);

    ps.growTo(2);
    ps[0] = Lit(b.ini(s, i));
    ps[1] = ~Lit(c.ini(s, i));
    s.addClause(ps);

    ps.growTo(3);
    ps[0] = ~Lit(a.ini(s, i));
    ps[1] = ~Lit(b.ini(s, i));
    ps[2] = Lit(c.ini(s, i));
    s.addClause(ps);
  }
}

void post_setunion(Solver &s, setvar a, setvar b, setvar c)
{
  for(int i = c.umin(s); i < std::min(a.umin(s), b.umin(s)); ++i)
    c.exclude(s, i, NO_REASON);
  for(int i = std::max(a.umax(s), b.umax(s))+1; i <= c.umax(s); ++i)
    c.exclude(s, i, NO_REASON);

  vec<Lit> ps;
  for(int i = c.umin(s); i <= c.umax(s); ++i) {
    ps.clear();

    if( i >= a.umin(s) && i <= a.umax(s) ) {
      ps.growTo(2);
      ps[0] = ~Lit(a.ini(s, i));
      ps[1] = Lit(c.ini(s, i));
      s.addClause(ps);
    }

    if( i >= b.umin(s) && i <= b.umax(s) ) {
      ps.growTo(2);
      ps[0] = ~Lit(b.ini(s, i));
      ps[1] = Lit(c.ini(s, i));
      s.addClause(ps);
    }

    ps.clear();
    pushifdef(ps, Lit(a.ini(s, i)));
    pushifdef(ps,  Lit(b.ini(s, i)));
    ps.push(~Lit(c.ini(s, i)));
    s.addClause(ps);
  }
}

void post_setsubseteq(Solver &s, setvar a, setvar b)
{
  for(int i = a.umin(s); i < b.umin(s); ++i)
    a.exclude(s, i, NO_REASON);
  for(int i = b.umax(s)+1; i <= a.umax(s); ++i)
    a.exclude(s, i, NO_REASON);

  vec<Lit> ps;
  for(int i = std::max(a.umin(s), b.umin(s)),
        iend = std::min(a.umax(s), b.umax(s))+1;
      i != iend; ++i) {
    ps.growTo(2);
    ps[0] = ~Lit( a.ini(s, i) );
    ps[1] = Lit( b.ini(s, i) );
    s.addClause(ps);
  }
}

void post_setsubset(Solver &s, setvar a, setvar b)
{
  post_setsubseteq(s, a, b);
  post_setneq(s, a, b);
}

void post_setsuperseteq(Solver &s, setvar a, setvar b)
{
  post_setsubseteq(s, b, a);
}

void post_setsuperset(Solver &s, setvar a, setvar b)
{
  post_setsubset(s, b, a);
}

void post_setsubseteq_re(Solver &s, setvar a, setvar b, Lit p)
{
  vec<Lit> ps;
  for(int i = a.umin(s); i < std::min(b.umin(s), a.umax(s)+1); ++i) {
    ps.growTo(2);
    ps[0] = ~Lit( a.ini(s, i) );
    ps[1] = ~p;
    s.addClause(ps);
  }
  for(int i = std::max(b.umax(s)+1, a.umin(s)); i <= a.umax(s); ++i) {
    ps.growTo(2);
    ps[0] = ~Lit( a.ini(s, i) );
    ps[1] = ~p;
    s.addClause(ps);
  }

  for(int i = std::max(a.umin(s), b.umin(s)),
        iend = std::min(a.umax(s), b.umax(s))+1;
      i < iend; ++i) {
    ps.clear();
    ps.push( ~Lit( a.ini(s, i) ) );
    pushifdef(ps, Lit( b.ini(s, i) ));
    ps.push(~p);
    s.addClause(ps);
  }

  ps.clear();
  vec<Lit> ps2;
  Var onlya = var_Undef;
  setneq::post_unique_used(s, a, b, ps, ps2, onlya);
  ps.clear();
  ps2.clear();
  if( onlya != var_Undef )
    ps2.push( Lit(onlya) );
  for(int i = std::max(a.umin(s), b.umin(s)),
        iend = std::min(a.umax(s), b.umax(s))+1;
      i < iend; ++i) {
    Var y = s.newVar();
    ps.clear();
    ps.push( Lit( a.ini(s, i) ) );
    ps.push( Lit(y) );
    s.addClause(ps);

    ps.clear();
    ps.push( ~Lit( b.ini(s, i) ) );
    ps.push( Lit(y) );
    s.addClause(ps);

    ps2.push( ~Lit(y) );
  }
  ps2.push(p);
  s.addClause(ps2);
}

void post_setsubseteq_re(Solver &s, setvar a, setvar b, cspvar r)
{
  assert( r.min(s) >= 0 && r.max(s) <= 1);
  if( r.max(s) == 0 )
    post_setsubseteq_re(s, a, b, ~Lit( r.eqi(s, 0) ));
  else
    post_setsubseteq_re(s, a, b, Lit( r.eqi(s, 1) ));
}

void post_setsuperseteq_re(Solver &s, setvar a, setvar b, Lit p)
{
  post_setsubseteq_re(s, b, a, p);
}


void post_setsuperseteq_re(Solver &s, setvar a, setvar b, cspvar r)
{
  post_setsubseteq_re(s, b, a, r);
}

} //namespace minicsp
