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

#include <algorithm>
#include "solver.hpp"
#include "cons.hpp"

using std::cout;
using std::vector;
using std::min;
using std::max;
using std::sort;
using std::ostream;

namespace minicsp {

/* This is based on Beldiceanu's algorithm from CP 2001, but only for
   pruning N. His algorithm for pruning X ends up being quadratic in
   the worst case (needs to examine every clique for every variable)
   and enforces RC. We do only BC here and we can reduce this to
   another sort+linear scan.
 */
class cons_atmostnvalue : public cons
{
  vector<cspvar> _x;
  vector<cspvar> _x2;
  cspvar _N;

  int dmin, dmax;
public:
  cons_atmostnvalue(Solver &s, vector<cspvar> const& x, cspvar N) :
    _x(x), _x2(x), _N(N)
  {
    dmin = x[0].min(s);
    dmax = x[0].max(s);
    for(size_t i = 1; i != _x.size(); ++i) {
      dmin = min( dmin, x[i].min(s) );
      dmax = max( dmax, x[i].max(s) );
      s.schedule_on_lb(_x[i], this);
      s.schedule_on_ub(_x[i], this);
    }
    s.schedule_on_ub(_N, this);

    DO_OR_THROW(propagate(s));
  }

  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
  Clause *propagate(Solver &s);
  void clone(Solver &other) {
    cons *con = new cons_atmostnvalue(other, _x, _N);
    other.addConstraint(con);
  }
};

namespace {
  struct incr_lb {
    Solver &s;
    incr_lb(Solver &ps) : s(ps) {}
    bool operator()(cspvar x1, cspvar x2) {
      return x1.min(s) < x2.min(s);
    }
  };

  struct decr_ub {
    Solver &s;
    decr_ub(Solver &ps) : s(ps) {}
    bool operator()(cspvar x1, cspvar x2) {
      return x1.max(s) > x2.max(s);
    }
  };
}

ostream& cons_atmostnvalue::print(Solver& s, ostream& os) const
{
  os << "atmostnvalue([";
  for(size_t i = 0; i != _x.size(); ++i) {
    if( i ) os << ", ";
    os << cspvar_printer(s, _x[i]);
  }
  os << "], " << cspvar_printer(s, _N) << ")";
  return os;
}

ostream& cons_atmostnvalue::printstate(Solver&s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != _x.size(); ++i) {
    os << cspvar_printer(s, _x[i]) << " in " << domain_as_range(s, _x[i])
       << ", ";
  }
  os << cspvar_printer(s, _N) << " in " << domain_as_range(s, _N);
  os << ")";
  return os;
}

/* The clauses generated here are as follows.

   --------------------------------------------------

   For pruning N:

   When we set N >= k, We only need to demontrate k distinct
   domains. We use the var that gave us the previous 'high' at the
   moment that reinit becomes true. This is guaranteed to have a
   greater lower bound than the variable that started this clique and
   has the smallest upper bound of the clique, so it is definitely
   suitable. We use both its lower and upper bound in the reason,
   unless it belongs to the first or last clique.

   --------------------------------------------------

   For pruning X:

   We get the set of literals produced above plus:

   - N <= ub(N)

   - for the first & last clique, we add the literals describing the
     lower (upper) bound of the highvar, which was excluded from ps
     when pruning N.

   - for every clique, we need not only that it is disconnected from
     the next, but also the particular lower and upper bounds of the
     common interval. We have a description of the upper bound
     (because we use highvar). If the min of high var is equal to the
     low of the clique, we do not add anything. Otherwise, we add
     *both* the min and max of the variable with the highest min in
     the clique. We need both, because we want to show membership in
     the clique and the existence of the min.

   Note that the literals of the lowvar are removed after we finish
   with each clique.

 */
Clause* cons_atmostnvalue::propagate(Solver &s)
{
  const size_t m = _x.size();

  sort(_x.begin(), _x.end(), incr_lb(s));

  int lb = 1;
  int reinit = 1;
  int low = dmin-1, high = dmax+1;

  vec<Lit> ps;
  vector<int> lows;
  vector<int> highs;
  vector<int> lowvars, highvars;

  lows.push_back( _x[0].min(s) );
  highs.push_back( _x[0].max(s) );
  lowvars.push_back( 0 );
  highvars.push_back( 0 );

  int highvar = 1;
  for(size_t i = 1; i < m; ) {
    bool newlow = false, newhigh = false;
    i = i + 1 - reinit;
    if( low < _x[i-1].min(s) ) {
      low = _x[i-1].min(s);
      newlow = true;
    }
    if( high > _x[i-1].max(s) ) {
      high = _x[i-1].max(s);
      highvar = i-1;
      newhigh = true;
    }
    reinit = (low > high);
    if( reinit ) {
      if( lows.size() > 1 )
        pushifdef( ps, _x[highvar].r_min(s) );
      ps.push( _x[highvar].r_max(s) );

      low = _x[i-1].min(s);
      high = _x[i-1].max(s);
      highvar = i-1;
      lows.push_back(low);
      highs.push_back(high);
      lowvars.push_back(i-1);
      highvars.push_back(i-1);
    } else {
      lows.back() = low;
      highs.back() = high;
      if( newlow )
        lowvars.back() = i-1;
      if( newhigh )
        highvars.back() = i-1;
    }
    lb += reinit;
  }
  pushifdef( ps, _x[highvar].r_min(s) );

  DO_OR_RETURN( _N.setminf(s, lb, ps) );

  if( lb < _N.max(s) )
    return 0L; // no other pruning possible

  pushifdef( ps, _N.r_max(s) );

  // note that none of the prunings below can cause a wipeout,
  // otherwise we would have detected failure above
  size_t interval = 0;
  int litsadded = 0;
  litsadded += pushifdef(ps, _x[highvars[0]].r_min(s));
  if( _x[ highvars[0]].min(s) < lows[0] ) {
    litsadded += pushifdef(ps, _x[lowvars[0]].r_min(s));
    litsadded += pushifdef(ps, _x[lowvars[0]].r_max(s));
  }
  for(size_t i = 0; i != m; ++i) {
    if( _x[i].min(s) > highs[interval] ) {
      if( _x[highvars[interval]].min(s) < lows[interval] ) {
        // we added something
        while( litsadded-- )
          ps.pop();
      }
      ++interval;
      litsadded = 0;
      if( _x[highvars[interval]].min(s) < lows[interval] ) {
        // the variable from that clique that is already in ps does
        // not explain the low of the interval
        litsadded += pushifdef(ps, _x[lowvars[interval]].r_min(s));
        litsadded += pushifdef(ps, _x[lowvars[interval]].r_max(s));
      } else {
        ;
      }
    }
    if( interval == 0 )
      _x[i].setminf(s, lows[interval], ps );
    else if(_x[i].min(s) > highs[interval-1] ) {
      bool tmpadded = pushifdef( ps, _x[i].r_geq(s, highs[interval-1]+1) );
      _x[i].setminf(s, lows[interval], ps );
      if( tmpadded )
        ps.pop();
    }
  }

  sort(_x2.begin(), _x2.end(), decr_ub(s));
  interval=0;
  size_t ni = highs.size();
  litsadded = 0;
  litsadded += pushifdef(ps, _x[highvars[ni-1]].r_max(s));
  if( _x[ highvars[ni-1]].min(s) < lows[0] ) {
    litsadded += pushifdef(ps, _x[lowvars[ni-1]].r_min(s));
    litsadded += pushifdef(ps, _x[lowvars[ni-1]].r_max(s));
  }
  for(size_t i = 0; i != m; ++i) {
    if( _x2[i].max(s) < lows[ni-1-interval] ) {
      int realiv = ni-1-interval;
      if( _x[highvars[realiv]].min(s) < lows[realiv] ) {
        // we added something
        while( litsadded-- )
          ps.pop();
      }
      ++interval;
      --realiv;
      if( _x[highvars[realiv]].min(s) < lows[realiv] ) {
        // the variable from that clique that is already in ps does
        // not explain the low of the interval
        litsadded = 0;
        litsadded += pushifdef(ps, _x[lowvars[realiv]].r_min(s));
        litsadded += pushifdef(ps, _x[lowvars[realiv]].r_max(s));
      }
    }
    if( interval == 0 )
      _x2[i].setmaxf(s, highs[ni-1-interval], ps );
    else if( _x2[i].max(s) < lows[ni-interval] ) {
      bool tmpadded = pushifdef( ps, _x2[i].r_leq(s, lows[ni-interval]-1) );
      _x2[i].setmaxf(s, highs[ni-1-interval], ps );
      if( tmpadded )
        ps.pop();
    }
  }

  return 0L;
}

void post_atmostnvalue(Solver &s, std::vector<cspvar> const& x, cspvar N)
{
  cons *con = new cons_atmostnvalue(s, x, N);
  s.addConstraint(con);
}

/* Independent set constraint: Given an adjacency matrix of a graph,
   compute a lower bound on the size of the maximum independent
   set.

   This uses the lower bound described in

   M. Halld\'orsson and J. Radhakrishnan. Greed is good: Approximating
   independent sets in sparse and bounded-degree graphs. In
   Proceedings STOC-94, pages 439{448, 1994.

   Roughly: remove the vertex of minimum degree along with its
   neighborhood. Repeat until no vertices remain. This is a lower
   bound on the size of the maximum independent set.

*/
class cons_independent_set : public cons
{
  int _n;
  vector<Var> _x;
  cspvar _N;

  size_t idx(size_t i, size_t j);
public:
  cons_independent_set(Solver &s, size_t n, std::vector<Var> const& x,
                       cspvar N) :
    _n(n), _x(x), _N(N)
  {
    for(size_t i = 0; i != _x.size(); ++i)
      s.schedule_on_lit(_x[i], this);
  }

  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
  Clause *propagate(Solver &s);
  void clone(Solver &other) {
    cons *con = new cons_independent_set(other, _n, _x, _N);
    other.addConstraint(con);
  }
};

ostream& cons_independent_set::print(Solver& s, ostream& os) const
{
  os << "independent_set([";
  for(size_t i = 0; i != _x.size(); ++i) {
    if( i ) os << ", ";
    os << var_printer(s, _x[i]);
  }
  os << "], " << cspvar_printer(s, _N) << ")";
  return os;
}

ostream& cons_independent_set::printstate(Solver&s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != _x.size(); ++i) {
    os << lit_printer(s, Lit(_x[i])) << ", ";
  }
  os << cspvar_printer(s, _N) << " in " << domain_as_range(s, _N);
  os << ")";
  return os;
}

size_t cons_independent_set::idx(size_t i, size_t j)
{
  if( i > j ) return idx(j, i);
  return _n*(_n-1)/2 - (_n-i)*(_n-i-1)/2 + j-i-1;
}

Clause *cons_independent_set::propagate(Solver &s)
{
  int nv = _n;
  int lb = 0;
  vec<Lit> ps;
  vector<int> degrees(_n, 0);
  vector<unsigned char> removed(_n, 0);
  while( nv > 0 ) {
    ++lb;
    for(int i = 0; i != _n; ++i) {
      if( removed[i] ) continue;
      for(int j = i+1; j != _n; ++j) {
        if( removed[j] ) continue;
        if( s.value( _x[idx(i,j)]) != l_False ) {
          ++degrees[i];
          ++degrees[j];
        }
      }
    }

    int t = -1; // toremove
    int md = degrees[0];
    for(int i = 0; i != _n; ++i) {
      if( removed[i] ) continue;
      if(t < 0 || md > degrees[i]) {
        t = i;
        md = degrees[i];
      }
    }

    assert(!removed[t]);
    removed[t] = 1;
    --nv;
    for(int i = 0; i != _n; ++i) {
      if( i == t ) continue;
      if( removed[i] ) continue;
      if( s.value( _x[idx(i,t)] ) != l_False ) {
        removed[i] = 1;
        --nv;
      } else {
        ps.push( Lit( _x[idx(i,t)] ) );
        degrees[i] = 0;
      }
    }
  }
  return _N.setminf(s, lb, ps);
}

void post_independent_set(Solver &s, size_t n, std::vector<Var> const& x,
                          cspvar N)
{
  cons *con = new cons_independent_set(s, n, x, N);
  s.addConstraint(con);
}

void post_atmostnvalue_md(Solver &s, std::vector<cspvar> const& x, cspvar N)
{
  vector<Var> eq_re;
  for(size_t i = 0; i != x.size(); ++i)
    for(size_t j = i+1; j != x.size(); ++j) {
      Var v = s.newVar();
      eq_re.push_back(v);
      post_eq_re(s, x[i], x[j], 0, Lit(v));
    }
  post_independent_set(s, x.size(), eq_re, N);
}

} //namespace minicsp
