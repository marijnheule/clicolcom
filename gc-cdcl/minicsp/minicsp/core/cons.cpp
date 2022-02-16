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
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <limits>
#include <cassert>
#include "cons.hpp"
#include "solver.hpp"

using std::ostream;
using std::pair;
using std::make_pair;
using std::vector;
using std::set;
using std::map;

// for temporary debugging
using std::cout;
#ifdef INVARIANTS
#define DOUT if(1) std::cout
#define DEXEC if(1)
#else
#define DOUT if(0) std::cout
#define DEXEC if(0)
#endif

namespace minicsp {

/* x == y + c     and        x == -y + c

   implemented as x = W*y + c, but with W == 1 or W == -1
*/
namespace eq {
  template<int W>
  Clause *eq_propagate(Solver &s, cspvar _x, cspvar _y, int _c, Lit event,
                       vec<Lit>& _reason)
  {
    domevent e = s.event(event);
    _reason.push(~event);
    switch(e.type) {
    case domevent::EQ:
      if( e.x == _x ) {
        return _y.assignf(s, W*e.d-W*_c, _reason);
      } else {
        return _x.assignf(s, W*e.d+_c, _reason);
      }
      break;
    case domevent::NEQ:
      if( e.x == _x ) {
        return _y.removef(s, W*e.d-W*_c, _reason);
      } else {
        return _x.removef(s, W*e.d+_c, _reason);
      }
      break;
    case domevent::GEQ:
      if( e.x == _x ) {
        if( W > 0 ) {
          return _y.setminf(s, e.d-_c, _reason);
        } else {
          return _y.setmaxf(s, -e.d+_c, _reason);
        }
      } else {
        if( W > 0 ) {
          return _x.setminf(s, e.d+_c, _reason);
        } else {
          return _x.setmaxf(s, -e.d+_c, _reason);
        }
      }
      break;
    case domevent::LEQ:
      if( e.x == _x ) {
        if( W > 0 ) {
          return _y.setmaxf(s, e.d-_c, _reason);
        } else {
          _reason.push(_y.e_geq(s, -e.d+_c));
          return _y.setminf(s, -e.d+_c, _reason);
        }
      } else {
        if( W > 0 ) {
          return _x.setmaxf(s, e.d+_c, _reason);
        } else {
          return _x.setminf(s, -e.d+_c, _reason);
        }
      }
      break;
    case domevent::NONE:
      assert(0);
    }
    return 0L;
  }

  template<int W>
  Clause *eq_initialize(Solver &s, cspvar _x, cspvar _y, int _c,
                        vec<Lit>& _reason)
  {
    if( W > 0 ) {
      {
        PUSH_TEMP( _reason,  _x.r_min(s) );
        DO_OR_RETURN(_y.setminf(s, _x.min(s)-_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _x.r_max(s) );
        DO_OR_RETURN(_y.setmaxf(s, _x.max(s)-_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _y.r_min(s) );
        DO_OR_RETURN(_x.setminf(s, _y.min(s)+_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _y.r_max(s) );
        DO_OR_RETURN(_x.setmaxf(s, _y.max(s)+_c, _reason));
      }
    } else {
      {
        PUSH_TEMP( _reason,  _x.r_max(s) );
        DO_OR_RETURN(_y.setminf(s, -_x.max(s)+_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _x.r_min(s) );
        DO_OR_RETURN(_y.setmaxf(s, -_x.min(s)+_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _y.r_max(s) );
        DO_OR_RETURN(_x.setminf(s, -_y.max(s)+_c, _reason));
      }

      {
        PUSH_TEMP( _reason,  _y.r_min(s) );
        DO_OR_RETURN(_x.setmaxf(s, -_y.min(s)+_c, _reason));
      }
    }

    for(int i = _x.min(s)+1, iend = _x.max(s); i < iend; ++i) {
      if( !_x.indomain(s, i) )
        DO_OR_RETURN(eq_propagate<W>(s, _x, _y, _c, _x.e_neq(s, i),
                                     _reason));
    }

    for(int i = _y.min(s)+1, iend = _y.max(s); i < iend; ++i) {
      if( !_y.indomain(s, i) )
        DO_OR_RETURN(eq_propagate<W>(s, _x, _y, _c, _y.e_neq(s, i),
                                     _reason));
    }

    return 0L;
  }
}

template<int W>
class cons_eq : public cons {
  cspvar _x, _y;
  int _c;
  vec<Lit> _reason; // cache to avoid the cost of allocation

  void requires() {
#ifndef _MSC_VER
    int w[ (W==1)+(W==-1)-1 ] __attribute__ ((__unused__));
#endif
  }
public:
  cons_eq(Solver &s,
           cspvar x, cspvar y, int c) :
    _x(x), _y(y), _c(c)
  {
    requires();
    // wake on everything! more immediate consequences
    s.wake_on_dom(x, this);
    s.wake_on_lb(x, this);
    s.wake_on_ub(x, this);
    s.wake_on_fix(x, this);
    s.wake_on_dom(y, this);
    s.wake_on_lb(y, this);
    s.wake_on_ub(y, this);
    s.wake_on_fix(y, this);
    _reason.capacity(2);

    if( eq::eq_initialize<W>(s, _x, _y, _c, _reason) )
      throw unsat();
  }

  virtual Clause *wake(Solver& s, Lit p) {
    _reason.clear();
    return eq::eq_propagate<W>(s, _x, _y, _c, p, _reason);
  }
  virtual void clone(Solver& other);
  virtual ostream& print(Solver &s, ostream& os) const;
  virtual ostream& printstate(Solver& s, ostream& os) const;
};

template<int W>
void cons_eq<W>::clone(Solver &other)
{
  cons *con = new cons_eq<W>(other, _x, _y, _c);
  other.addConstraint(con);
}

template<int W>
ostream& cons_eq<W>::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " = ";
  if( W == -1 ) os << "-";
  os << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  return os;
}

template<int W>
ostream& cons_eq<W>::printstate(Solver& s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  os << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x) << ", "
     << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ")";
  return os;
}

/* x == y + c */
void post_eq(Solver& s, cspvar x, cspvar y, int c)
{
  cons *con = new cons_eq<1>(s, x, y, c);
  s.addConstraint(con);
}

/* x == -y + c */
void post_neg(Solver &s, cspvar x, cspvar y, int c)
{
  cons *con = new cons_eq<-1>(s, x, y, c);
  s.addConstraint(con);
}

/* x != y + c */
namespace neq {
  Clause *neq_propagate(Solver &s, cspvar _x, cspvar _y, int _c, Lit event,
                        vec<Lit>& _reason);

  Clause *neq_initialize(Solver &s, cspvar _x, cspvar _y, int _c,
                       vec<Lit>& _reason)
  {
    if( _x.min(s) == _x.max(s) )
      return neq_propagate(s, _x, _y, _c,
                           _x.e_eq(s, _x.min(s)),
                           _reason);
    else if( _y.min(s) == _y.max(s) )
      return neq_propagate(s, _x, _y, _c,
                           _y.e_eq(s, _y.min(s)),
                           _reason);
    return 0L;
  }

  Clause *neq_propagate(Solver &s, cspvar _x, cspvar _y, int _c, Lit event,
                       vec<Lit>& _reason)
  {
    domevent e = s.event(event);
    if( e.type != domevent::EQ )
      return 0L;
    _reason.push(~event);
    if( e.x == _x ) {
      return _y.removef(s, e.d-_c, _reason);
    } else {
      return _x.removef(s, e.d+_c, _reason);
    }
    return 0L;
  }
}

class cons_neq : public cons {
  cspvar _x, _y;
  int _c;
  vec<Lit> _reason; // cache to avoid the cost of allocation
public:
  cons_neq(Solver &s,
           cspvar x, cspvar y, int c) :
    _x(x), _y(y), _c(c)
  {
    s.wake_on_fix(x, this);
    s.wake_on_fix(y, this);
    _reason.capacity(2);
    DO_OR_THROW(neq::neq_initialize(s, _x, _y, _c, _reason));
  }

  virtual Clause *wake(Solver& s, Lit p) override;
  virtual void clone(Solver& other) override;
  virtual ostream& print(Solver &s, ostream& os) const override;
  virtual ostream &printstate(Solver &s, ostream &os) const override;
};

Clause *cons_neq::wake(Solver& s, Lit event)
{
  _reason.clear();
  return neq::neq_propagate(s, _x, _y, _c, event, _reason);
}

void cons_neq::clone(Solver &other)
{
  cons *con = new cons_neq(other, _x, _y, _c);
  other.addConstraint(con);
}

ostream& cons_neq::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " != " << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  return os;
}

ostream& cons_neq::printstate(Solver &s, ostream& os) const
{
  print(s, os);
  os << " (with " << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x)
     << ", " << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ")";
  return os;
}

/* x != y + c */
void post_neq(Solver& s, cspvar x, cspvar y, int c)
{
  cons *con = new cons_neq(s, x, y, c);
  s.addConstraint(con);
}

/* x = y + c <=> b */
namespace eq_re {
  /* Return true if the domains of _x and _y are disjoint (offset
     _c). Also describe the prunings that made it so in _reason. Note
     if it returns false, _reason contains garbage. */
  bool disjoint_domains(Solver &s, cspvar _x, cspvar _y, int _c,
                        int& watch,
                        vec<Lit>& _reason) {
    using std::min;
    using std::max;
    if( _x.min(s) > _y.max(s) + _c ) {
      pushifdef(_reason, _x.r_min(s));
      pushifdef(_reason, _y.r_max(s));
      return true;
    } else if( _x.max(s) < _y.min(s) + _c ) {
      pushifdef(_reason, _x.r_max(s));
      pushifdef(_reason, _y.r_min(s));
      return true;
    } else {
      if( _x.indomain(s, watch) && _y.indomain(s, watch - _c) )
        return false;
      if( _x.min(s) > _y.min(s) + _c )
        pushifdef(_reason, _x.r_min(s) );
      else
        pushifdef(_reason, _y.r_min(s) );
      if( _x.max(s) < _y.max(s) + _c )
        pushifdef(_reason, _x.r_max(s));
      else
        pushifdef(_reason, _y.r_max(s));
      for(int i = max(_x.min(s), _y.min(s) + _c),
            iend = min(_x.max(s), _y.max(s) + _c)+1; i != iend; ++i) {
        if( _x.indomain(s, i) ) {
          if( _y.indomain(s, i-_c) ) {
            watch = i;
            return false;
          } else
            pushifdef(_reason, _y.r_neq(s, i-_c) );
        } else
          pushifdef(_reason, _x.r_neq(s,i));
      }
      return true;
    }
  }
}

class cons_eq_re : public cons {
  cspvar _x, _y;
  int _c;
  Lit _b;
  vec<Lit> _reason;
  int _watch; // last known common value
public:
  cons_eq_re(Solver &s,
             cspvar x, cspvar y, int c,
             Lit b,
             int initwatch) :
    _x(x), _y(y), _c(c), _b(b),
    _watch(initwatch)
  {
    s.wake_on_dom(_x, this);
    s.wake_on_fix(_x, this);
    s.wake_on_dom(_y, this);
    s.wake_on_fix(_y, this);
    s.wake_on_lit(var(_b), this);
    _reason.capacity(5);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause *cons_eq_re::wake(Solver &s, Lit p)
{
  domevent pevent = s.event(p);
  _reason.clear();
  if( s.value(_b) == l_True ) {
    _reason.push( ~_b );
    if( p == _b )
      return eq::eq_initialize<1>(s, _x, _y, _c, _reason);
    else
      return eq::eq_propagate<1>(s, _x, _y, _c, p, _reason);
  } else if( s.value(_b) == l_False ) {
    _reason.push( _b );
    if( p == ~_b )
      return neq::neq_initialize(s, _x, _y, _c, _reason);
    else
      return neq::neq_propagate(s, _x, _y, _c, p, _reason);
  }

  if( ((pevent.type == domevent::NEQ && pevent.d == _watch) ||
       (pevent.type == domevent::GEQ && pevent.d > _watch) ||
       (pevent.type == domevent::LEQ && pevent.d < _watch)) &&
      eq_re::disjoint_domains(s, _x, _y, _c, _watch, _reason) ) {
    DO_OR_RETURN(s.enqueueFill(~_b, _reason));
  } else if( _x.min(s) == _y.min(s)+_c &&
             _x.min(s) == _x.max(s) &&
             _x.min(s) == _y.max(s)+_c ) {
    _reason.clear(); // disjoint_domains() may have put garbage here
    _reason.push(_x.r_eq(s));
    _reason.push(_y.r_eq(s));
    DO_OR_RETURN(s.enqueueFill(_b, _reason));
  }

  return 0L;
}

void cons_eq_re::clone(Solver & other)
{
  cons *con = new cons_eq_re(other, _x, _y, _c, _b, _watch);
  other.addConstraint(con);
}

ostream& cons_eq_re::print(Solver &s, ostream& os) const
{
  os << "(" << cspvar_printer(s, _x) << " = " << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  os << ") <=> " << lit_printer(s, _b);
  return os;
}

ostream& cons_eq_re::printstate(Solver & s, ostream& os) const
{
  print(s, os);
  os << " (with " << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x)
     << ", " << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ", " << lit_printer(s, _b)
     << ")";
  return os;
}

void post_eqneq_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  vec<Lit> reason;
  int dummywatch = x.min(s)-1;
  if( eq_re::disjoint_domains(s, x, y, c, dummywatch, reason) ) {
    s.uncheckedEnqueue(~b); // will not cause unsat, because b is
                            // unset by if test in callers
    return;
  } else if( x.min(s) == y.min(s)+c &&
             x.min(s) == x.max(s) &&
             x.min(s) == y.max(s)+c ) {
    s.uncheckedEnqueue(b);
    return;
  }

  cons *con = new cons_eq_re(s, x, y, c, b, dummywatch);
  s.addConstraint(con);
}

void post_eq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1 );
  if( b.min(s) == 1 ) {
    post_eq(s, x, y, c);
    return;
  } else if( b.max(s) == 0 ) {
    post_neq(s, x, y, c);
    return;
  }

  post_eqneq_re(s, x, y, c, Lit(b.eqi(s, 1)));
}

void post_eq_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  if( s.value(b) == l_True ) {
    post_eq(s, x, y, c);
    return;
  } else if( s.value(b) == l_False ) {
    post_neq(s, x, y, c);
    return;
  }
  post_eqneq_re(s, x, y, c, b);
}

void post_neq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1 );
  if( b.min(s) == 1 ) {
    post_neq(s, x, y, c);
    return;
  } else if( b.max(s) == 0 ) {
    post_eq(s, x, y, c);
    return;
  }

  post_eqneq_re(s, x, y, c, ~Lit(b.eqi(s, 1)));
}

void post_neq_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  if( s.value(b) == l_True ) {
    post_neq(s, x, y, c);
    return;
  } else if( s.value(b) == l_False ) {
    post_eq(s, x, y, c);
    return;
  }
  post_eqneq_re(s, x, y, c, ~b);
}

/* x <= y + c */
namespace leq {
  Clause *leq_propagate(Solver &s, cspvar _x, cspvar _y, int _c,
                        vec<Lit>& _reason)
  {
    if( _y.max(s) + _c < _x.min(s) ) { // failure
      pushifdef(_reason, _x.r_min(s));
      pushifdef(_reason, _y.r_leq( s, _x.min(s) - _c - 1));
      Clause *r = Clause_new(_reason);
      s.addInactiveClause(r);
      return r;
    }

    if( _y.min(s) + _c < _x.min(s) ) {
      pushifdef(_reason, _x.r_min(s));
      _y.setminf(s, _x.min(s) - _c, _reason);
    }
    if( _x.max(s) - _c > _y.max(s) ) {
      pushifdef(_reason, _y.r_max(s));
      _x.setmaxf(s, _y.max(s) + _c, _reason);
    }

    return 0L;
  }
}

class cons_le : public cons {
  cspvar _x, _y;
  int _c;
  vec<Lit> _reason; // cache to avoid the cost of allocation
public:
  cons_le(Solver &s,
          cspvar x, cspvar y, int c) :
    _x(x), _y(y), _c(c)
  {
    s.wake_on_lb(x, this);
    s.wake_on_ub(y, this);
    _reason.push(); _reason.push();
    if( wake(s, lit_Undef) )
      throw unsat();
  }

  virtual Clause *wake(Solver& s, Lit p);
  virtual void clone(Solver& othersolver);
  virtual ostream& print(Solver &s, ostream& os) const;
};

Clause *cons_le::wake(Solver& s, Lit)
{
  _reason.clear();
  return leq::leq_propagate(s, _x, _y, _c, _reason);
}

void cons_le::clone(Solver &other)
{
  cons *con = new cons_le(other, _x, _y, _c);
  other.addConstraint(con);
}

ostream& cons_le::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " <= " << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  return os;
}

// v1 <= v2 + c
void post_leq(Solver& s, cspvar v1, cspvar v2, int c)
{
  cons *con = new cons_le(s, v1, v2, c);
  s.addConstraint(con);
}

// v1 < v2 + c
void post_less(Solver& s, cspvar v1, cspvar v2, int c)
{
  cons *con = new cons_le(s, v1, v2, c-1);
  s.addConstraint(con);
}

/* (x <= y + c) <=> b and the half-implication cases */
template<bool LI, bool RI>
class cons_leq_re : public cons {
  cspvar _x, _y;
  int _c;
  Lit _b;
  vec<Lit> _reason;
public:
  cons_leq_re(Solver &s,
              cspvar x, cspvar y, int c, Lit b) :
    _x(x), _y(y), _c(c), _b(b)
  {
    s.wake_on_lb(_x, this);
    s.wake_on_ub(_x, this);
    s.wake_on_lb(_y, this);
    s.wake_on_ub(_y, this);
    s.wake_on_lit(var(_b), this);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

template<bool LI, bool RI>
Clause *cons_leq_re<LI, RI>::wake(Solver &s, Lit)
{
  _reason.clear();

  if( LI && s.value(_b) == l_True ) {
    _reason.push( ~_b );
    return leq::leq_propagate(s, _x, _y, _c, _reason);
  } else if( RI && s.value(_b) == l_False ) {
    _reason.push( _b );
    return leq::leq_propagate(s, _y, _x, -_c-1, _reason);
  }

  if( LI && _x.min(s) > _y.max(s) + _c ) {
    pushifdef(_reason, _x.r_min(s));
    pushifdef(_reason, _y.r_max(s));
    DO_OR_RETURN(s.enqueueFill(~_b, _reason));
  } else if( RI && _x.max(s) <= _y.min(s) + _c ) {
    pushifdef(_reason, _x.r_max(s));
    pushifdef(_reason, _y.r_min(s));
    DO_OR_RETURN(s.enqueueFill(_b, _reason));
  }

  return 0L;
}

template<bool LI, bool RI>
void cons_leq_re<LI, RI>::clone(Solver & other)
{
  cons *con = new cons_leq_re<LI, RI>(other, _x, _y, _c, _b);
  other.addConstraint(con);
}

template<bool LI, bool RI>
ostream& cons_leq_re<LI, RI>::print(Solver &s, ostream& os) const
{
  os << "(" << cspvar_printer(s, _x) << " <= " << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  os << ") ";
  if( LI ) os << '<';
  os << '=';
  if( RI ) os << '>';
  os << lit_printer(s, _b);
  return os;
}

template<bool LI, bool RI>
ostream& cons_leq_re<LI, RI>::printstate(Solver & s, ostream& os) const
{
  print(s, os);
  os << " (with " << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x)
     << ", " << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ", " << lit_printer(s, _b)
     << ")";
  return os;
}

template<bool LI, bool RI>
void post_leq_re_common(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  if( LI && s.value(b) == l_True ) {
    post_leq(s, x, y, c);
    return;
  } else if( RI && s.value(b) == l_False ) {
    post_less(s, y, x, -c);
    return;
  }

  if( RI && x.max(s) <= y.min(s) + c ) {
    s.uncheckedEnqueue(b);
    return;
  } else if( LI && x.min(s) > y.max(s) + c ) {
    s.uncheckedEnqueue(~b);
    return;
  }

  cons *con = new cons_leq_re<LI, RI>(s, x, y, c, b);
  s.addConstraint(con);
}

void post_leq_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
    post_leq_re_common<true, true>(s, x, y, c, b);
}

void post_leq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1 );
  if( b.max(s) == 0 ) {
    post_leq_re(s, x, y, c, ~Lit(b.eqi(s, 0)));
    return;
  }

  post_leq_re(s, x, y, c, Lit(b.eqi(s, 1)));
}

void post_less_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re(s, x, y, c-1, b);
}

void post_less_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re(s, x, y, c-1, b);
}

void post_geq_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re(s, y, x, -c, b);
}

void post_geq_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re(s, y, x, -c, b);
}

void post_gt_re(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re(s, y, x, -c-1, b);
}

void post_gt_re(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re(s, y, x, -c-1, b);
}

void post_leq_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
    post_leq_re_common<false, true>(s, x, y, c, b);
}

void post_leq_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1 );
  if( b.max(s) == 0 ) {
    post_leq_re_ri(s, x, y, c, ~Lit(b.eqi(s, 0)));
    return;
  }

  post_leq_re_ri(s, x, y, c, Lit(b.eqi(s, 1)));
}

void post_less_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_ri(s, x, y, c-1, b);
}

void post_less_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_ri(s, x, y, c-1, b);
}

void post_geq_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_ri(s, y, x, -c, b);
}

void post_geq_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_ri(s, y, x, -c, b);
}

void post_gt_re_ri(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_ri(s, y, x, -c-1, b);
}

void post_gt_re_ri(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_ri(s, y, x, -c-1, b);
}

void post_leq_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
    post_leq_re_common<true, false>(s, x, y, c, b);
}

void post_leq_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  assert( b.min(s) >= 0 && b.max(s) <= 1 );
  if( b.max(s) == 0 ) {
    post_leq_re_li(s, x, y, c, ~Lit(b.eqi(s, 0)));
    return;
  }

  post_leq_re_li(s, x, y, c, Lit(b.eqi(s, 1)));
}

void post_less_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_li(s, x, y, c-1, b);
}

void post_less_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_li(s, x, y, c-1, b);
}

void post_geq_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_li(s, y, x, -c, b);
}

void post_geq_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_li(s, y, x, -c, b);
}

void post_gt_re_li(Solver &s, cspvar x, cspvar y, int c, cspvar b)
{
  post_leq_re_li(s, y, x, -c-1, b);
}

void post_gt_re_li(Solver &s, cspvar x, cspvar y, int c, Lit b)
{
  post_leq_re_li(s, y, x, -c-1, b);
}

/* abs: |x| = y + c*/
class cons_abs : public cons {
  cspvar _x, _y;
  int _c;
  vec<Lit> _reason; // cache to avoid the cost of allocation
public:
  cons_abs(Solver &s,
           cspvar x, cspvar y, int c) :
    _x(x), _y(y), _c(c)
  {
    using std::max;
    using std::min;

    // we only wake on dom here, otherwise it is too complicated
    s.wake_on_dom(x, this);
    s.wake_on_dom(y, this);
    _reason.capacity(3); _reason.growTo(2, lit_Undef);

    if( y.min(s) + _c < 0 )
      DO_OR_THROW(_y.setmin(s, 0, NO_REASON));

    if( _y.min(s) + _c > max(abs(x.max(s)), abs(x.min(s))) ||
        (x.min(s) > 0 &&
         _y.max(s) + _c < min(abs(x.max(s)), abs(x.min(s)))) )
      throw unsat();

    for(int i = _x.min(s), iend = _x.max(s)+1; i != iend; ++i) {
      if( !_x.indomain(s, i) )
        DO_OR_THROW(wake(s, _x.e_neq(s, i)));
      if( !_y.indomain(s, abs(i)-_c) )
        _x.remove(s, i, NO_REASON);
    }

    for(int i = _y.min(s), iend = _y.max(s)+1; i != iend; ++i) {
      if( !_y.indomain(s, i) )
        DO_OR_THROW(wake(s, _y.e_neq(s, i)));
      if( !_x.indomain(s, i+_c) && !_x.indomain(s, -i-_c) )
        _y.remove(s, i, NO_REASON);
    }
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause* cons_abs::wake(Solver &s, Lit event)
{
  domevent e = s.event(event);
  _reason[0] = ~event;

  if( e.x == _x ) {
    if( _y.indomain(s, abs(e.d)-_c) && !_x.indomain(s, -e.d) ) {
      _reason[1] = _x.r_neq( s, -e.d );
      if( var(_reason[1]) == var_Undef ) {
        _reason[1] = _y.e_neq(s, abs(e.d)-_c);
        return _y.remove(s, abs(e.d)-_c, _reason);
      }
      _reason.push( _y.e_neq( s, abs(e.d)-_c) );
      Clause *confl = _y.remove(s, abs(e.d)-_c, _reason);
      _reason.shrink(1);
      return confl;
    }
  } else {
    if( e.d < 0 ) return 0L;
    if( _x.indomain(s, e.d+_c) ) {
      _reason[1] = _x.e_neq( s, e.d+_c);
      DO_OR_RETURN(_x.remove(s, e.d+_c, _reason));
    }
    if( _x.indomain(s, -e.d-_c) ) {
      _reason[1] = _x.e_neq( s, -e.d-_c );
      DO_OR_RETURN(_x.remove(s, -e.d-_c, _reason));
    }
  }

  return 0L;
}

void cons_abs::clone(Solver &other)
{
  cons *con = new cons_abs(other, _x, _y, _c);
  other.addConstraint(con);
}

ostream& cons_abs::print(Solver &s, ostream& os) const
{
  os << "abs(" << cspvar_printer(s, _x) << ") = "
     << cspvar_printer(s, _y);
  if( _c > 0 )
    os << " + " << _c;
  else if( _c < 0 )
    os << " - " << -_c;
  return os;
}

ostream& cons_abs::printstate(Solver& s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  os << cspvar_printer(s, _x) << " in " << domain_as_range(s, _x) << ", "
     << cspvar_printer(s, _y) << " in " << domain_as_range(s, _y) << ")";
  return os;
}

void post_abs(Solver& s, cspvar v1, cspvar v2, int c)
{
  if( v1.min(s) >= 0 ) {
    post_eq(s, v1, v2, c);
    return;
  }
  if( v1.max(s) <= 0 ) {
    post_neg(s, v2, v1, -c);
    return;
  }
  cons *con = new cons_abs(s, v1, v2, c);
  s.addConstraint(con);
}

/* cons_lin_le

   N-ary linear inequality

   sum w[i]*v[i] + c <= 0

   The template parameter N is used to tell the compiler to create
   optimized versions for small N. When N = 0, it uses the generic
   loop, otherwise constant propagation will make it a fixed loop.
*/
template<size_t N>
class cons_lin_le : public cons {
  int _c;

  vec<Lit> _ps;
  size_t n;
  vector<int> pspos;
  pair<int, cspvar> _vars[0];

  // computes the minimum contribution of var x with coeff w:
  // w*x.min() if w>0, w*x.max() if w<0
  int vmin(Solver &s, int w, cspvar x);
  // describe the domain of x to use in the reason for whatever we do:
  // 'x>=x.min()' if w>0, 'x<=x.max()' if w<0
  Lit litreason(Solver &s, int lb, int w, cspvar x);
public:
  cons_lin_le(Solver &s,
              vector< pair<int, cspvar> > const& vars,
              int c);

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;

  void dispose() { this->~cons_lin_le(); free(this); }
};

template<size_t N>
cons *new_cons_lin_le(Solver &s,
                      vector< pair<int, cspvar> > const& vars,
                      int c)
{
  void *mem = malloc(sizeof(cons_lin_le<N>) +
                     vars.size()*sizeof( pair<int, cspvar> ));
  return new (mem) cons_lin_le<N>(s, vars, c);
}

template<size_t N>
int cons_lin_le<N>::vmin(Solver &s, int w, cspvar x)
{
  if( w > 0 ) return w*x.min(s);
  else return w*x.max(s);
}

template<size_t N>
Lit cons_lin_le<N>::litreason(Solver &s, int lb,
                              int w, cspvar x)
{
  if( w > 0 )
    return x.r_geq(s, x.min(s));
  else
    return x.r_leq(s, x.max(s));
}

template<size_t N>
cons_lin_le<N>::cons_lin_le(Solver &s,
                            vector< pair<int, cspvar> > const& vars,
                            int c) :
  _c(c)
{
  assert(N == 0 || N == vars.size());
  n = vars.size();
  int lb = _c;
  for(size_t i = 0; i != vars.size(); ++i) {
    _vars[i] = vars[i];
    lb += vmin(s, _vars[i].first, _vars[i].second);
    if( _vars[i].first > 0 )
      s.wake_on_lb(_vars[i].second, this);
    else
      s.wake_on_ub(_vars[i].second, this);
  }
  if( lb > 0 ) throw unsat();
  _ps.growTo(n, lit_Undef);
  if( wake(s, lit_Undef) )
    throw unsat();
}

template<size_t N>
Clause *cons_lin_le<N>::wake(Solver &s, Lit)
{
  pspos.clear();
  pspos.resize(n);

  int lb = _c;
  size_t nl = 0;
  for(size_t i = 0; i != n; ++i) {
    int w = _vars[i].first;
    cspvar x = _vars[i].second;
    _ps[nl] = litreason(s, lb, w, x);
    pspos[i] = nl;
    if( toInt(_ps[nl]) >= 0 ) ++nl;
    lb += vmin(s, w, x);
  }

  if( lb > 0 ) {
    // shrink _ps to contruct the clause, then bring it back to its
    // previous size
    _ps.shrink(n - nl );
    Clause *r = Clause_new(_ps);
    _ps.growTo(n, lit_Undef);
    s.addInactiveClause(r);
    return r;
  }

  _ps.shrink(n - nl );
  for(size_t i = 0; i != n; ++i) {
    Lit l = _ps[pspos[i]];
    int w = _vars[i].first;
    cspvar x = _vars[i].second;
    int gap = lb-vmin(s, w, x);
    if( w > 0 ) {
      int ibound = -gap/w;
      if( litreason(s, lb, w, x) == lit_Undef ) {
        _ps.push(x.e_leq(s, ibound));
        x.setmax(s, ibound, _ps);
        _ps.pop();
      } else {
        _ps[pspos[i]] = x.e_leq(s, ibound);
        x.setmax(s, ibound, _ps);
        _ps[pspos[i]] = l;
      }
    } else {
      int ibound;
      // rounding towards zero is weird
      if( gap < 0 )
        ibound = -gap/w;
      else
        ibound = -(gap-w-1)/w;
      if( litreason(s, lb, w, x) == lit_Undef ) {
        _ps.push(x.e_geq(s, ibound));
        x.setmin(s, ibound, _ps);
        _ps.pop();
      } else {
        _ps[pspos[i]] = x.e_geq(s, ibound);
        x.setmin(s, ibound, _ps);
        _ps[pspos[i]] = l;
      }
    }
  }
  _ps.growTo(n, lit_Undef);

  return 0L;
}

template<size_t N>
void cons_lin_le<N>::clone(Solver &other)
{
  vector< pair<int, cspvar> > v(_vars, _vars+n);
  cons *con = new_cons_lin_le<N>(other, v, _c);
  other.addConstraint(con);
}

template<size_t N>
ostream& cons_lin_le<N>::print(Solver &s, ostream& os) const
{
  for(size_t i = 0; i != n; ++i) {
    if( _vars[i].first == 1 ) {
      if( i != 0 )
        os << " + ";
      os << cspvar_printer(s, _vars[i].second);
    } else if( _vars[i].first == -1 ) {
      if( i != 0 )
        os << " ";
      os << "- " << cspvar_printer(s, _vars[i].second);
    } else if( _vars[i].first > 0 ) {
      if( i != 0 )
        os << " +";
      os << _vars[i].first << "*" << cspvar_printer(s, _vars[i].second);
    }
    else if( _vars[i].first < 0 ) {
      if( i != 0 )
        os << " ";
      os << "- " << -_vars[i].first
         << "*" << cspvar_printer(s, _vars[i].second);
    }
  }
  if( _c > 0 )
    os << " + " << _c;
  else
    os << " - " << -_c;
  os << " <= 0";

  return os;
}

template<size_t N>
ostream& cons_lin_le<N>::printstate(Solver& s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != n; ++i) {
    if( i ) os << ", ";
    cspvar x = _vars[i].second;
    os << cspvar_printer(s, x) << " in " << domain_as_range(s, x);
  }
  os << ")";
  return os;
}

namespace lin {
  struct cmp_varid {
    bool operator()(pair<int, cspvar> x1, pair<int, cspvar> x2) const {
      return x1.second.id() < x2.second.id();
    }
  };
}

void post_lin_leq(Solver &s, vector<cspvar> const& vars,
                   vector<int> const &coeff, int c)
{
  assert(vars.size() == coeff.size());
  vector< pair<int, cspvar> > pairs;
  for(size_t i = 0; i != vars.size(); ++i) {
    if( !coeff[i] ) continue;
    if( vars[i].min(s) == vars[i].max(s) ) {
      c += coeff[i]*vars[i].min(s);
      continue;
    }
    pairs.push_back( make_pair(coeff[i], vars[i]));
  }

  if( pairs.size() == 0 ) {
    if( c > 0 ) throw unsat();
    return;
  }

  sort(pairs.begin(), pairs.end(), lin::cmp_varid());
  vector< pair<int, cspvar> >::iterator i = pairs.begin(), j = i;
  ++j;
  for(; j != pairs.end(); ++j) {
    if( i->second == j->second )
      i->first += j->first;
    else {
      ++i;
      *i = *j;
    }
  }
  ++i;
  pairs.erase(i, pairs.end());

  cons *con;
  switch(pairs.size()) {
  case 1:   con = new_cons_lin_le<1>(s, pairs, c); break;
  case 2:   con = new_cons_lin_le<2>(s, pairs, c); break;
  case 3:   con = new_cons_lin_le<3>(s, pairs, c); break;
  default:  con = new_cons_lin_le<0>(s, pairs, c); break;
  }
  s.addConstraint(con);
}

void post_lin_less(Solver &s, vector<cspvar> const& vars,
                    vector<int> const &coeff, int c)
{
  post_lin_leq(s, vars, coeff, c+1);
}

void post_lin_leq_right_imp_re(Solver &s, vector<cspvar> const&vars,
                               vector<int> const &coeff,
                               int c, cspvar b)
{
  assert(vars.size() == coeff.size());
  assert(b.min(s) >= 0 && b.max(s) <= 1);

  if(b.max(s) == 0) {
    vector<int> c1(coeff);
    for(size_t i = 0; i != vars.size(); ++i)
      c1[i] = -c1[i];
    post_lin_leq(s, vars, c1, -c+1);
    return;
  }

  vector<cspvar> v1(vars);
  vector<int> c1(coeff);

  /* Let L = c + sum coeff[i]*vars[i] */
  /* Discover lb, ub, s.t. lb <= L <= ub */
  int ub = c, lb = c;
  for(size_t i = 0; i != vars.size(); ++i) {
    if( coeff[i] > 0 ) {
      ub += coeff[i] * vars[i].max(s);
      lb += coeff[i] * vars[i].min(s);
    } else {
      ub += coeff[i] * vars[i].min(s);
      lb += coeff[i] * vars[i].max(s);
    }
    c1[i] = -c1[i];
  }

  if( ub <= 0 ) { /* L <= 0 always, b is true */
    Clause *confl = b.assign(s, 1, NO_REASON);
    if( confl ) throw unsat();
    return;
  }
  if( lb > 0 ) { /* L > 0 always, rhs unaffected */
    return;
  }

  v1.push_back(b);
  c1.push_back(lb-1);

  // post -L + 1 + (lb - 1)*b <= 0
  post_lin_leq(s, v1, c1, -c+1);
}

void post_lin_less_right_imp_re(Solver &s, vector<cspvar> const&vars,
                                vector<int> const &coeff,
                                int c, cspvar b)
{
  post_lin_leq_right_imp_re(s, vars, coeff, c+1, b);
}

void post_lin_leq_left_imp_re(Solver &s,
                              std::vector<cspvar> const&vars,
                              std::vector<int> const &coeff,
                              int c,
                              cspvar b)
{
  assert(vars.size() == coeff.size());
  assert(b.min(s) >= 0 && b.max(s) <= 1);

  if( b.min(s) == 1 ) {
    post_lin_leq(s, vars, coeff, c);
    return;
  }

  vector<cspvar> v1(vars);
  vector<int> c1(coeff);

  /* Let L = c + sum coeff[i]*vars[i] */
  /* Discover lb, ub, s.t. lb <= L <= ub */
  int ub = c, lb = c;
  for(size_t i = 0; i != vars.size(); ++i) {
    if( coeff[i] > 0 ) {
      ub += coeff[i] * vars[i].max(s);
      lb += coeff[i] * vars[i].min(s);
    } else {
      ub += coeff[i] * vars[i].min(s);
      lb += coeff[i] * vars[i].max(s);
    }
  }

  if( ub <= 0 ) { /* L <= 0 always, lhs unaffected */
    return;
  }
  if( lb > 0 ) { /* L > 0 always, set b to false */
    Clause *confl = b.assign(s, 0, NO_REASON);
    if( confl ) throw unsat();
    return;
  }

  v1.push_back(b);
  c1.push_back(ub+1);

  // post -L + 1 + (lb - 1)*b <= 0
  post_lin_leq(s, v1, c1, c-ub-1);
}

void post_lin_less_left_imp_re(Solver &s,
                               std::vector<cspvar> const&vars,
                               std::vector<int> const &coeff,
                               int c,
                               cspvar b)
{
  post_lin_leq_left_imp_re(s, vars, coeff, c+1, b);
}

void post_lin_leq_iff_re(Solver &s, std::vector<cspvar> const& vars,
                          std::vector<int> const& coeff,
                          int c, cspvar b)
{
  post_lin_leq_left_imp_re(s, vars, coeff, c, b);
  post_lin_leq_right_imp_re(s, vars, coeff, c, b);
}

void post_lin_less_iff_re(Solver &s, std::vector<cspvar> const& vars,
                           std::vector<int> const& coeff,
                           int c, cspvar b)
{
  post_lin_leq_left_imp_re(s, vars, coeff, c, b);
  post_lin_leq_right_imp_re(s, vars, coeff, c, b);
}

// sum coeff[i]*vars[i] + c = 0
void post_lin_eq(Solver &s,
                 std::vector<cspvar> const& vars,
                 std::vector<int> const& coeff,
                 int c)
{
  post_lin_leq(s, vars, coeff, c);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_leq(s, vars, c1, -c);
}

// sum coeff[i]*vars[i] + c = 0 implies b = 1
void post_lin_eq_right_imp_re(Solver &s,
                              std::vector<cspvar> const& vars,
                              std::vector<int> const& coeff,
                              int c,
                              cspvar b)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_leq_right_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_leq_right_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps;
  ps.push( ~Lit(b1.eqi(s, 1)) );
  ps.push( ~Lit(b2.eqi(s, 1)) );
  ps.push( Lit(b.eqi(s, 1)) );
  s.addClause(ps);
}

// b = 1 implies sum coeff[i]*vars[i] + c = 0
void post_lin_eq_left_imp_re(Solver &s,
                             std::vector<cspvar> const& vars,
                             std::vector<int> const& coeff,
                             int c,
                             cspvar b)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_leq_left_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_leq_left_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps;
  ps.push( ~Lit(b1.eqi(s, 1)) );
  ps.push( ~Lit(b2.eqi(s, 1)) );
  ps.push( Lit(b.eqi(s, 1)) );
  s.addClause(ps);
}

// b=1 iff sum coeff[i]*vars[i] + c = 0
//
// L >= 0 => b1, L <= 0 => b2, b <=> b1 /\ b2
void post_lin_eq_iff_re(Solver &s,
                        std::vector<cspvar> const& vars,
                        std::vector<int> const& coeff,
                        int c,
                        cspvar b)
{
  assert(b.min(s)>=0 && b.max(s) <=1);
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_leq_iff_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_leq_iff_re(s, vars, c1, -c, b2);

  Lit bl;
  if( b.max(s) == 0 )
    bl = ~Lit(b.eqi(s, 0));
  else
    bl = Lit(b.eqi(s, 1));

  vec<Lit> ps1, ps2, ps3;
  ps1.push( ~bl );
  ps1.push( Lit( b1.eqi(s, 1) ) );

  ps2.push( ~bl );
  ps2.push( Lit( b2.eqi(s, 1) ) );

  ps3.push( ~Lit( b1.eqi(s, 1) ) );
  ps3.push( ~Lit( b2.eqi(s, 1) ) );
  ps3.push( bl );

  s.addClause(ps1);
  s.addClause(ps2);
  s.addClause(ps3);
}

/* linear inequality: L != 0

   implemented as
   L >= 0 => b1, L <= 0 => b2, not b1 or not b2 */
void post_lin_neq(Solver &s,
                  std::vector<cspvar> const& vars,
                  std::vector<int> const& coeff,
                  int c)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_leq_right_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_leq_right_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps;
  ps.push( ~Lit(b1.eqi(s, 1)) );
  ps.push( ~Lit(b2.eqi(s, 1)) );
  s.addClause(ps);
}

/* L != 0 => b = 1

   implemented as
   L > 0 => b1, L < 0 => b2, b1 or b2 => b
*/
void post_lin_neq_right_imp_re(Solver &s,
                               std::vector<cspvar> const& vars,
                               std::vector<int> const& coeff,
                               int c,
                               cspvar b)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_less_right_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_less_right_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps;
  ps.push( ~Lit(b1.eqi(s, 1)) );
  ps.push( ~Lit(b2.eqi(s, 1)) );
  ps.push( Lit( b.eqi(s, 1)) );
  s.addClause(ps);
}

/* b = 1 => L != 0

   implemented as
   b1 => L > 0, b2 => L < 0, b => b1 or b2
*/
void post_lin_neq_left_imp_re(Solver &s,
                              std::vector<cspvar> const& vars,
                              std::vector<int> const& coeff,
                              int c,
                              cspvar b)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_less_left_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_less_left_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps;
  ps.push( Lit(b1.eqi(s, 1)) );
  ps.push( Lit(b2.eqi(s, 1)) );
  ps.push( ~Lit( b.eqi(s, 1)) );
  s.addClause(ps);
}

/* b = 1 <=> L != 0

   implemented as
   b1 => L > 0, b2 => L < 0, b <=> b1 or b2
 */
void post_lin_neq_iff_re(Solver &s,
                         std::vector<cspvar> const& vars,
                         std::vector<int> const& coeff,
                         int c,
                         cspvar b)
{
  cspvar b1 = s.newCSPVar(0,1);
  cspvar b2 = s.newCSPVar(0,1);
  post_lin_less_left_imp_re(s, vars, coeff, c, b1);
  vector<int> c1(coeff);
  for(size_t i = 0; i != vars.size(); ++i)
    c1[i] = -coeff[i];
  post_lin_less_left_imp_re(s, vars, c1, -c, b2);
  vec<Lit> ps1, ps2, ps3;
  ps1.push( Lit(b1.eqi(s, 1)) );
  ps1.push( Lit(b2.eqi(s, 1)) );
  ps1.push( ~Lit( b.eqi(s, 1)) );

  ps2.push( ~Lit(b2.eqi(s, 1)) );
  ps2.push( Lit( b.eqi(s, 1)) );

  ps3.push( ~Lit(b2.eqi(s, 1)) );
  ps3.push( Lit( b.eqi(s, 1)) );

  s.addClause(ps1);
  s.addClause(ps2);
  s.addClause(ps3);
}

/* cons_pb

   PseudoBoolean constraint

   sum_i weight[i]*var[i] >= lb
*/
namespace pb {
  struct compare_abs_weights {
    bool operator()(pair<int, Var> a, pair<int, Var> b)
    {
      return abs(a.first) > abs(b.first);
    }
  };
}

class cons_pb : public cons {
  std::vector< std::pair<int, Var> > _posvars;
  std::vector< std::pair<int, Var> > _negvars;
  std::vector< std::pair<int, Var> > _svars; // sorted by absolute value
  int _lb;
  int _ub; // the perfect upper bound
public:
  cons_pb(Solver &s,
          std::vector<Var> const& vars,
          std::vector<int> const& weights, int lb);

  virtual Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

void post_pb(Solver& s, vector<Var> const& vars,
              vector<Var> const& weights, int lb)
{
  cons *c = new cons_pb(s, vars, weights, lb);
  s.addConstraint(c);
}

void post_pb(Solver& s, vector<cspvar> const& vars,
              vector<Var> const& weights, int lb)
{
  vector<Var> vbool(vars.size());
  for(size_t i = 0; i != vars.size(); ++i)
    vbool[i] = vars[i].eqi(s, 1);

  cons *c = new cons_pb(s, vbool, weights, lb);
  s.addConstraint(c);
}

// sum ci*bi >= lb ==> b
void post_pb_right_imp_re(Solver& s, std::vector<cspvar> const& vars,
                          std::vector<int> const& weights, int lb,
                          cspvar b)
{
  vector<cspvar> v1(vars);
  vector<int> w1(weights);
  int lub = 0, llb = 0;
  for(size_t i = 0; i != vars.size(); ++i) {
    if( weights[i] > 0 )
      lub += weights[i];
    else
      llb += weights[i];
    w1[i] = -w1[i];
  }

  int diff = lub - llb + 1;
  v1.push_back(b);
  w1.push_back( diff );
  post_pb(s, v1, w1, - lb + 1);
}

// b ==> sum ci*bi >= lb
void post_pb_left_imp_re(Solver& s, std::vector<cspvar> const& vars,
                         std::vector<int> const& weights, int lb,
                         cspvar b)
{
  vector<cspvar> v1(vars);
  vector<int> w1(weights);
  int lub = 0, llb = 0;
  for(size_t i = 0; i != vars.size(); ++i) {
    if( weights[i] > 0 )
      lub += weights[i];
    else
      llb += weights[i];
  }

  int diff = lub - llb + 1;
  v1.push_back(b);
  w1.push_back( -diff );
  post_pb(s, v1, w1, lb - diff);
}

void post_pb_iff_re(Solver& s, std::vector<cspvar> const& vars,
                    std::vector<int> const& weights, int lb,
                    cspvar b)
{
  post_pb_right_imp_re(s, vars, weights, lb, b);
  post_pb_left_imp_re(s, vars, weights, lb, b);
}

cons_pb::cons_pb(Solver& s,
                 vector<Var> const& vars, vector<int> const& weights,
                 int lb) : _lb(lb), _ub(0)
{
  assert(vars.size() == weights.size());
  size_t n = vars.size();
  for(size_t i = 0; i != n; ++i) {
    if( weights[i] > 0 ) {
      _posvars.push_back(make_pair(weights[i], vars[i]));
      _ub += weights[i];
    }
    else if( weights[i] < 0 )
      _negvars.push_back(make_pair(weights[i], vars[i]));
    else if( weights[i] == 0 ) continue;
    _svars.push_back(make_pair(weights[i], vars[i]));
  }
  sort(_svars.begin(), _svars.end(), pb::compare_abs_weights());

  // wake on every assignment
  n = _svars.size();
  for(size_t i = 0; i != n; ++i)
    s.wake_on_lit(_svars[i].second, this);
}

Clause *cons_pb::wake(Solver& s, Lit)
{
  size_t np = _posvars.size(),
    nn = _negvars.size();
  int ub = 0;
  for(size_t i = 0; i != np; ++i) {
    pair<int, Var> const& cv = _posvars[i];
    int q = toInt(s.value(cv.second));
    /* the expression 1 - ((-q+1) >> 1) is
       1 <=> true or unset
       0 <=> false */
    int p = 1-((1-q)>>1);
    ub += p*cv.first;
  }
  for(size_t i = 0; i != nn; ++i) {
    pair<int, Var> const& cv = _negvars[i];
    int q = toInt(s.value(cv.second));
    /* the expression (q+1)>>1 is
       0 <=> false or unset
       1 <=> true */
    int p = (q+1)>>1;
    ub += p*cv.first;
  }
  if( ub >= _lb ) return 0L; // satisfiable

  vec<Lit> ps;
  size_t n = _svars.size();
  int aub = _ub;
  for(size_t i = 0; i != n; ++i) {
    pair<int, Var> const& cv = _svars[i];
    lbool sv = s.value(cv.second);
    if( sv == l_Undef ) continue;
    if( cv.first > 0 && sv == l_False ) {
      aub -= cv.first;
      Lit lx(cv.second);
      ps.push( lx );
    } else if( cv.first < 0 && sv == l_True ) {
      aub += cv.first;
      Lit lx(cv.second, true);
      ps.push( lx );
    }
    if( aub < _lb ) break;
  }

  Clause *conf = Clause_new(ps, true);
  s.addInactiveClause(conf);
  return conf;
}

void cons_pb::clone(Solver& other)
{
  vector<Var> v;
  vector<int> w;
  for(size_t i = 0; i != _svars.size(); ++i) {
    v.push_back(_svars[i].second);
    w.push_back(_svars[i].first);
  }
  cons *con = new cons_pb(other, v, w, _lb);
  other.addConstraint(con);
}

ostream& cons_pb::print(Solver &s, ostream& os) const
{
  for(size_t i = 0; i != _svars.size(); ++i) {
    if( _svars[i].first == 1 ) {
      if( i != 0 )
        os << " + ";
      os << var_printer(s, _svars[i].second);
    } else if( _svars[i].first == -1 ) {
      if( i != 0 )
        os << " ";
      os << "- " << var_printer(s, _svars[i].second);
    } else if( _svars[i].first > 0 ) {
      if( i != 0 )
        os << " +";
      os << _svars[i].first << "*" << var_printer(s, _svars[i].second);
    }
    else if( _svars[i].first < 0 ) {
      if( i != 0 )
        os << " ";
      os << "- " << -_svars[i].first << "*"
         << var_printer(s, _svars[i].second);
    }
  }
  os << " >= " << _lb;
  return os;
}

ostream& cons_pb::printstate(Solver& s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != _svars.size(); ++i) {
    if( i ) os << ", ";
    Var x = _svars[i].second;
    int xmin = 0, xmax = 1;
    if( s.value(x) == l_True ) xmin = 1;
    else if( s.value(x) == l_False ) xmax = 0;
    os << var_printer(s, x) << " in [" << xmin << ", " << xmax << "]";
  }
  os << ")";
  return os;
}

/* pseudoboolean with a var on the rhs: sum w[i]*v[i] = rhs */
class cons_pbvar : public cons {
  std::vector< std::pair<int, Var> > _posvars;
  std::vector< std::pair<int, Var> > _negvars;
  std::vector< std::pair<int, Var> > _svars; // sorted by absolute value
  int _c;
  cspvar _rhs;
  vec<Lit> _ubreason, _lbreason;
public:
  cons_pbvar(Solver &s,
             std::vector< std::pair<int, Var> > const &vars,
             int c, cspvar rhs)
    : _c(c), _rhs(rhs)
  {
    size_t n = vars.size();
    int lb = _c, ub = _c;
    for(size_t i = 0; i != n; ++i) {
      if( vars[i].first > 0 ) {
        _posvars.push_back( vars[i] );
        ub += vars[i].first;
      }
      else if( vars[i].first < 0 ) {
        _negvars.push_back( vars[i] );
        lb += vars[i].first;
      }
      else
        continue;
      _svars.push_back( vars[i] );
    }
    sort(_svars.begin(), _svars.end(), pb::compare_abs_weights());

    _rhs.setmin(s, lb, NO_REASON);
    _rhs.setmax(s, ub, NO_REASON);

    // wake on every assignment
    n = _svars.size();
    for(size_t i = 0; i != n; ++i)
      s.wake_on_lit(_svars[i].second, this);
    s.wake_on_lb(_rhs, this);
    s.wake_on_ub(_rhs, this);

    _ubreason.growTo(n+1);
    _lbreason.growTo(n+1);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause *cons_pbvar::wake(Solver &s, Lit)
{
  size_t np = _posvars.size(),
    nn = _negvars.size(),
    n = np+nn;

  int lb = _c, ub = _c;

  _ubreason.growTo(n+1);
  _lbreason.growTo(n+1);

  size_t ubi = 0, lbi = 0;
  for(size_t i = 0; i != np; ++i) {
    pair<int, Var> const& cv = _posvars[i];
    int q = toInt(s.value(cv.second));
    if( q >= 0 )
      ub += cv.first;
    else
      _ubreason[ubi++] = Lit(cv.second);
    if( q > 0 ) {
      lb += cv.first;
      _lbreason[lbi++] = Lit(cv.second, true);
    }
  }
  for(size_t i = 0; i != nn; ++i) {
    pair<int, Var> const& cv = _negvars[i];
    int q = toInt(s.value(cv.second));
    if( q >= 0 )
      lb += cv.first;
    else
      _lbreason[lbi++] = Lit(cv.second);
    if( q > 0 ) {
      ub += cv.first;
      _ubreason[ubi++] = Lit(cv.second, true);
    }
  }

  bool nonimpliedlb = false,
    nonimpliedub = false;
  if( _rhs.min(s) > lb )
    nonimpliedlb = true;
  if( _rhs.max(s) < ub )
    nonimpliedub = true;

  _lbreason.shrink(_lbreason.size() - lbi);
  DO_OR_RETURN(_rhs.setminf(s, lb, _lbreason));
  _ubreason.shrink(_ubreason.size() - ubi);
  DO_OR_RETURN(_rhs.setmaxf(s, ub, _ubreason));

  // check again, because rhs might have had holes!
  if( _rhs.min(s) > lb )
    nonimpliedlb = true;
  if( _rhs.max(s) < ub )
    nonimpliedub = true;

  for(int i = 0; i != _ubreason.size(); ++i)
    _lbreason.push(_ubreason[i]);

  int rhslb = _rhs.min(s);
  int rhsub = _rhs.max(s);

  // note we gather everything in _lbreason now
  if( nonimpliedlb )
    _lbreason.push( _rhs.r_min(s) );
  if( nonimpliedub )
    _lbreason.push( _rhs.r_max(s) );

  for(size_t i = 0; i != np; ++i) {
    pair<int, Var> const& cv = _posvars[i];
    int q = toInt(s.value(cv.second));
    if( q ) continue;
    if( lb + cv.first > rhsub )
      DO_OR_RETURN(s.enqueueFill( Lit( cv.second, true ), _lbreason ));
    if( ub - cv.first < rhslb )
      DO_OR_RETURN(s.enqueueFill( Lit( cv.second ), _lbreason ));
  }
  for(size_t i = 0; i != nn; ++i) {
    pair<int, Var> const& cv = _negvars[i];
    int q = toInt(s.value(cv.second));
    if( q ) continue;
    if( ub + cv.first < rhslb )
      DO_OR_RETURN(s.enqueueFill( Lit(cv.second, true ), _lbreason ));
    if( lb - cv.first > rhsub )
      DO_OR_RETURN(s.enqueueFill( Lit(cv.second ), _lbreason ));
  }

  return 0L;
}

void cons_pbvar::clone(Solver &other)
{
  cons *con = new cons_pbvar(other, _svars, _c, _rhs);
  other.addConstraint(con);
}

ostream& cons_pbvar::print(Solver &s, ostream& os) const
{
  for(size_t i = 0; i != _svars.size(); ++i) {
    if( _svars[i].first == 1 ) {
      if( i != 0 )
        os << " + ";
      os << var_printer(s, _svars[i].second);
    } else if( _svars[i].first == -1 ) {
      if( i != 0 )
        os << " ";
      os << "- " << var_printer(s, _svars[i].second);
    } else if( _svars[i].first > 0 ) {
      if( i != 0 )
        os << " +";
      os << _svars[i].first << "*" << var_printer(s, _svars[i].second);
    }
    else if( _svars[i].first < 0 ) {
      if( i != 0 )
        os << " ";
      os << "- " << -_svars[i].first
         << "*" << var_printer(s, _svars[i].second);
    }
  }
  os << " = " << cspvar_printer(s, _rhs);

  return os;
}

ostream& cons_pbvar::printstate(Solver& s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != _svars.size(); ++i) {
    if( i ) os << ", ";
    Var x = _svars[i].second;
    int xmin = 0, xmax = 1;
    if( s.value(x) == l_True ) xmin = 1;
    else if( s.value(x) == l_False ) xmax = 0;
    os << var_printer(s, x) << " in [" << xmin << ", " << xmax << "]";
  }
  os << ", " << cspvar_printer(s, _rhs)
     << " in " << domain_as_range(s, _rhs);
  os << ")";
  return os;
}

void post_pb(Solver &s, std::vector<Var> const& vars,
             std::vector<int> const& weights, int c, cspvar rhs)
{
  vector< pair< int, Var> > vbool;
  vbool.reserve(vars.size());
  for(size_t i = 0; i != vars.size(); ++i) {
    if( s.value(vars[i]) == l_True ) {
      c += weights[i];
      continue;
    } else if( s.value(vars[i]) == l_False ) {
      continue;
    }
    vbool.push_back( make_pair( weights[i], vars[i] ));
  }

  cons *con = new cons_pbvar(s, vbool, c, rhs);
  s.addConstraint(con);
}

void post_pb(Solver& s, std::vector<cspvar> const& vars,
             std::vector<int> const& weights, int c, cspvar rhs)
{
  vector< pair< int, Var> > vbool;
  vbool.reserve(vars.size());
  for(size_t i = 0; i != vars.size(); ++i) {
    if( vars[i].min(s) == 1 ) {
      c += weights[i];
      continue;
    } else if( vars[i].max(s) == 0 ) {
      continue;
    }
    vbool.push_back( make_pair( weights[i], vars[i].eqi(s, 1)));
  }

  cons *con = new cons_pbvar(s, vbool, c, rhs);
  s.addConstraint(con);
}

/* x = y*z */
class cons_mult : public cons {
  cspvar _x, _y, _z;
  vec<Lit> _reason;

  // divide and round up
  int divup(int d, int n) {
    if( d > 0 ) {
      if( n > 0 )
        return (d+n-1)/n;
      else
        return d/n;
    } else {
      if( n > 0 )
        return d/n;
      else
        return (d-n-1)/n;
    }
  }
  // divide and round down
  int divdn(int d, int n) {
    if( d > 0 ) {
      if( n > 0 )
        return d/n;
      else
        return (d-n-1)/n;
    } else {
      if( n > 0 )
        return (d-n+1)/n;
      else
        return d/n;
    }
  }
public:
  cons_mult(Solver &s, cspvar x, cspvar y, cspvar z) :
    _x(x), _y(y), _z(z)
  {
    s.wake_on_lb(_x, this);
    s.wake_on_ub(_x, this);
    s.wake_on_lb(_y, this);
    s.wake_on_ub(_y, this);
    s.wake_on_lb(_z, this);
    s.wake_on_ub(_z, this);

    _reason.capacity(5);
    wake(s, lit_Undef);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause* cons_mult::wake(Solver &s, Lit)
{
  using std::min;
  using std::max;
  if( _y.min(s) > 0 && _z.min(s) > 0 ) {
    _reason.clear();
    pushifdef(_reason, _z.r_min(s));
    pushifdef(_reason, _y.r_min(s));
    DO_OR_RETURN(_x.setminf( s, _y.min(s)*_z.min(s), _reason ));

    _reason.clear();
    pushifdef(_reason, _z.r_max(s));
    pushifdef(_reason, _y.r_max(s));
    DO_OR_RETURN(_x.setmaxf( s, _y.max(s)*_z.max(s), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_min(s));
    pushifdef(_reason, _z.r_max(s));
    DO_OR_RETURN(_y.setminf( s, divup(_x.min(s), _z.max(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_max(s));
    pushifdef(_reason, _z.r_min(s));
    DO_OR_RETURN(_y.setmaxf( s, divdn(_x.max(s), _z.min(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_min(s));
    pushifdef(_reason, _y.r_max(s));
    DO_OR_RETURN(_z.setminf( s, divup(_x.min(s), _y.max(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_max(s));
    pushifdef(_reason, _y.r_min(s));
    DO_OR_RETURN(_z.setmaxf( s, divdn(_x.max(s), _y.min(s)), _reason ));
  } if( _y.max(s) < 0 && _z.max(s) < 0 ) {
    _reason.clear();
    pushifdef(_reason, _z.r_max(s));
    pushifdef(_reason, _y.r_max(s));
    DO_OR_RETURN(_x.setminf( s, _y.max(s)*_z.max(s), _reason ));

    _reason.clear();
    pushifdef(_reason, _z.r_min(s));
    pushifdef(_reason, _y.r_min(s));
    DO_OR_RETURN(_x.setmaxf( s, _y.min(s)*_z.min(s), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_max(s));
    pushifdef(_reason, _z.r_max(s));
    DO_OR_RETURN(_y.setminf( s, divup(_x.max(s), _z.max(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_min(s));
    pushifdef(_reason, _z.r_min(s));
    DO_OR_RETURN(_y.setmaxf( s, divdn(_x.min(s), _z.min(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_max(s));
    pushifdef(_reason, _y.r_max(s));
    DO_OR_RETURN(_z.setminf( s, divup(_x.max(s), _y.max(s)), _reason ));

    _reason.clear();
    pushifdef(_reason, _x.r_min(s));
    pushifdef(_reason, _y.r_min(s));
    DO_OR_RETURN(_z.setmaxf( s, divdn(_x.min(s), _y.min(s)), _reason ));
  } else {
    _reason.clear();
    pushifdef(_reason, _y.r_min(s));
    pushifdef(_reason, _y.r_max(s));
    pushifdef(_reason, _z.r_min(s));
    pushifdef(_reason, _z.r_max(s));
    DO_OR_RETURN(_x.setminf( s, min( min( _y.min(s)*_z.min(s),
                                          _y.min(s)*_z.max(s)),
                                     min( _y.max(s)*_z.min(s),
                                          _y.max(s)*_z.max(s))),
                             _reason));
    DO_OR_RETURN(_x.setmaxf( s, max( max( _y.min(s)*_z.min(s),
                                          _y.min(s)*_z.max(s)),
                                     max( _y.max(s)*_z.min(s),
                                          _y.max(s)*_z.max(s))),
                             _reason));
    if( _y.min(s) > 0 ) {
      _reason.clear();
      pushifdef(_reason, _y.r_min(s));
      pushifdef(_reason, _y.r_max(s));
      pushifdef(_reason, _x.r_min(s));
      pushifdef(_reason, _x.r_max(s));
      DO_OR_RETURN(_z.setminf( s, min( min( _x.min(s) / _y.min(s),
                                            _x.min(s) / _y.max(s)),
                                       min( _x.max(s) / _y.min(s),
                                            _x.max(s) / _y.max(s))),
                               _reason));
      DO_OR_RETURN(_z.setmaxf( s, max( max( _x.min(s) / _y.min(s),
                                            _x.min(s) / _y.max(s)),
                                       max( _x.max(s) / _y.min(s),
                                            _x.max(s) / _y.max(s))),
                               _reason));
    }
    if( _z.min(s) > 0 ) {
      _reason.clear();
      pushifdef(_reason, _z.r_min(s));
      pushifdef(_reason, _z.r_max(s));
      pushifdef(_reason, _x.r_min(s));
      pushifdef(_reason, _x.r_max(s));
      DO_OR_RETURN(_y.setminf( s, min( min( _x.min(s) / _z.min(s),
                                            _x.min(s) / _z.max(s)),
                                       min( _x.max(s) / _z.min(s),
                                            _x.max(s) / _z.max(s))),
                               _reason));
      DO_OR_RETURN(_y.setmaxf( s, max( max( _x.min(s) / _z.min(s),
                                            _x.min(s) / _z.max(s)),
                                       max( _x.max(s) / _z.min(s),
                                            _x.max(s) / _z.max(s))),
                               _reason));
    }
  }
  return 0L;
}

void cons_mult::clone(Solver& other)
{
  cons *con = new cons_mult(other, _x, _y, _z);
  other.addConstraint(con);
}

ostream& cons_mult::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " = "
     << cspvar_printer(s, _y) << "*" << cspvar_printer(s, _z);
  return os;
}

ostream& cons_mult::printstate(Solver &s, ostream& os) const
{
  print(s, os);
  os << " (with " <<  cspvar_printer(s, _x)
     << " in " << domain_as_range(s, _x)
     << ", " <<   cspvar_printer(s, _y)
     << " in " << domain_as_range(s, _y)
     << ", " <<  cspvar_printer(s, _z)
     << " in " << domain_as_range(s, _z);
  return os;
}

void post_mult(Solver& s, cspvar x, cspvar y, cspvar z)
{
  cons *con = new cons_mult(s, x, y, z);
  s.addConstraint(con);
}

void post_div(Solver& s, cspvar x, cspvar y, cspvar z)
{
  post_mult(s, y, x, z);
}

class cons_mod;

// x = min(y,z)
class cons_min : public cons {
  cspvar _x, _y, _z;
  int _c;
  vec<Lit> _reason;
public:
  cons_min(Solver &s,
           cspvar x, cspvar y, cspvar z) :
    _x(x), _y(y), _z(z)
  {
    s.wake_on_lb(_x, this);
    s.wake_on_ub(_x, this);
    s.wake_on_lb(_y, this);
    s.wake_on_ub(_y, this);
    s.wake_on_lb(_z, this);
    s.wake_on_ub(_z, this);
    wake(s, lit_Undef);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause *cons_min::wake(Solver &s, Lit)
{
  using std::min;

  _reason.clear();

  {
    PUSH_TEMP(_reason, _y.r_min(s));
    PUSH_TEMP(_reason, _z.r_min(s));
    DO_OR_RETURN(_x.setminf(s, min(_y.min(s), _z.min(s)), _reason));
  }
  {
    PUSH_TEMP(_reason, _y.r_max(s));
    PUSH_TEMP(_reason, _z.r_max(s));
    DO_OR_RETURN(_x.setmaxf(s, min(_y.max(s), _z.max(s)), _reason));
  }

  _reason.clear();
  pushifdef(_reason,_x.r_min(s));
  DO_OR_RETURN(_y.setminf(s, _x.min(s), _reason));
  DO_OR_RETURN(_z.setminf(s, _x.min(s), _reason));

  return 0L;
}

void cons_min::clone(Solver & other)
{
  cons *con = new cons_min(other, _x, _y, _z);
  other.addConstraint(con);
}

ostream& cons_min::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " = min("
     << cspvar_printer(s, _y) << ", "
     << cspvar_printer(s, _z) << ")";
  return os;
}

ostream& cons_min::printstate(Solver & s, ostream& os) const
{
  print(s, os);
  os << " (with " << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x)
     << ", " << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ", " << cspvar_printer(s, _z) << " in " << domain_as_set(s, _z)
     << ")";
  return os;
}

void post_min(Solver &s, cspvar x, cspvar y, cspvar z)
{
  cons *con = new cons_min(s, x, y, z);
  s.addConstraint(con);
}

// x = max(y,z)
class cons_max : public cons {
  cspvar _x, _y, _z;
  vec<Lit> _reason;
public:
  cons_max(Solver &s,
           cspvar x, cspvar y, cspvar z) :
    _x(x), _y(y), _z(z)
  {
    s.wake_on_lb(_x, this);
    s.wake_on_ub(_x, this);
    s.wake_on_lb(_y, this);
    s.wake_on_ub(_y, this);
    s.wake_on_lb(_z, this);
    s.wake_on_ub(_z, this);
    wake(s, lit_Undef);
  }

  Clause *wake(Solver& s, Lit p);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

Clause *cons_max::wake(Solver &s, Lit)
{
  using std::max;

  _reason.clear();

  {
    PUSH_TEMP(_reason, _y.r_min(s));
    PUSH_TEMP(_reason, _z.r_min(s));
    DO_OR_RETURN(_x.setminf(s, max(_y.min(s), _z.min(s)), _reason));
  }
  {
    PUSH_TEMP(_reason, _y.r_max(s));
    PUSH_TEMP(_reason, _z.r_max(s));
    DO_OR_RETURN(_x.setmaxf(s, max(_y.max(s), _z.max(s)), _reason));
  }

  _reason.clear();
  pushifdef(_reason,_x.r_max(s));
  DO_OR_RETURN(_y.setmaxf(s, _x.max(s), _reason));
  DO_OR_RETURN(_z.setmaxf(s, _x.max(s), _reason));

  return 0L;
}

void cons_max::clone(Solver & other)
{
  cons *con = new cons_max(other, _x, _y, _z);
  other.addConstraint(con);
}

ostream& cons_max::print(Solver &s, ostream& os) const
{
  os << cspvar_printer(s, _x) << " = max("
     << cspvar_printer(s, _y) << ", "
     << cspvar_printer(s, _z) << ")";
  return os;
}

ostream& cons_max::printstate(Solver & s, ostream& os) const
{
  print(s, os);
  os << " (with " << cspvar_printer(s, _x) << " in " << domain_as_set(s, _x)
     << ", " << cspvar_printer(s, _y) << " in " << domain_as_set(s, _y)
     << ", " << cspvar_printer(s, _z) << " in " << domain_as_set(s, _z)
     << ")";
  return os;
}

void post_max(Solver &s, cspvar x, cspvar y, cspvar z)
{
  cons *con = new cons_max(s, x, y, z);
  s.addConstraint(con);
}


/* Element: R = X[I]

   Note that indexing is 1-based here, following flatzinc's element */
namespace element {
  size_t idx(int i, int imin,
             int j, int rmin,
             int n) {
    return (j-rmin)*n + (i - imin);
  }
}

void post_element(Solver &s, cspvar R, cspvar I,
                  vector<cspvar> const& X,
                  int offset)
{
  using std::min;
  using std::max;

  assert(!X.empty());

  I.setmin(s, offset, NO_REASON);
  I.setmax(s, X.size()+offset-1, NO_REASON);

  int imin = I.min(s),
    imax = I.max(s);

  int rmin = X[imin-offset].min(s), rmax = X[imin-offset].max(s);
  for(int i = imin+1, iend = imax+1; i < iend; ++i) {
    rmin = min(rmin, X[i-offset].min(s));
    rmax = max(rmax, X[i-offset].max(s));
  }
  rmin = max(rmin, R.min(s));
  rmax = min(rmax, R.max(s));
  R.setmin(s, rmin, NO_REASON);
  R.setmax(s, rmax, NO_REASON);

  // -I=i -Xi=j R=j
  vec<Lit> ps;
  for(int i = imin; i <= imax; ++i) {
    for(int j = rmin; j <= rmax; ++j) {
      ps.clear();
      Lit xij = Lit(X[i-offset].eqi(s, j));
      if( xij == lit_Undef || s.value(xij) == l_False)
        continue; // true clause
      Lit rj = Lit(R.eqi(s, j));
      if( s.value(rj) == l_True )
        continue; // true clause
      ps.push( ~xij );
      ps.push( ~Lit(I.eqi(s, i)) );
      pushifdef( ps, rj );
      s.addClause(ps);
    }
  }

  int n = imax-imin+1;
  int d = rmax-rmin+1;

  Var yfirst, wfirst;

  yfirst = s.newVar();
  for(int i = 1; i != n*d; ++i)
    s.newVar();

  wfirst = s.newVar();
  for(int i = 1; i != n*d; ++i)
    s.newVar();

  using element::idx;

  vec<Lit> ps1, ps2, ps3;
  // define y_ij <=> R=j or Xi=j
  for(int i = imin; i <= imax; ++i)
    for(int j = rmin; j <= rmax; ++j) {
      ps1.clear(); ps2.clear(); ps3.clear();
      Var yij = yfirst + idx(i, imin, j, rmin, n);
      // R=j y_ij
      pushifdef(ps1, Lit(R.eqi(s, j)));
      ps1.push( Lit(yij) );
      s.addClause(ps1);

      // Xi=j y_ij
      pushifdef(ps2, Lit(X[i-offset].eqi(s, j)));
      ps2.push( Lit(yij) );
      s.addClause(ps2);

      // -y_ij -R=j -Xi=j
      Var rij, xij;
      rij = R.eqi(s, j);
      xij = X[i-offset].eqi(s, j);
      if( rij == var_Undef || xij == var_Undef )
        continue;
      ps3.push( ~Lit(yij) );
      ps3.push( ~Lit(rij) );
      ps3.push( ~Lit(xij) );
      s.addClause(ps3);
    }

  // -y_i1 ... -y_id -I=i
  for(int i = imin; i <= imax; ++i) {
    ps.clear();
    if( !I.indomain(s, i) ) continue;
    ps.push( ~Lit(I.eqi(s, i)));
    for(int j = rmin; j <= rmax; ++j) {
      Var yij = yfirst + idx(i, imin, j, rmin, n);
      ps.push( ~Lit(yij) );
    }
    s.addClause(ps);
  }

  // define w_ij <=> -I=i or -Xi=j
  for(int i = imin; i <= imax; ++i)
    for(int j = rmin; j <= rmax; ++j) {
      Var wij = wfirst + idx(i, imin, j, rmin, n);
      ps1.clear(); ps2.clear(); ps3.clear();
      // I=i, w_ij
      pushifdef( ps1, Lit(I.eqi(s,i)));
      ps1.push( Lit(wij) );
      s.addClause(ps1);

      // Xi=j, w_ij
      pushifdef( ps2, Lit(X[i-offset].eqi(s, j)));
      ps2.push( Lit(wij) );
      s.addClause(ps2);

      // -w_ij, -I=i, -Xi=j
      if( I.indomain(s, i) && X[i-offset].indomain(s, j) ) {
        ps3.push( ~Lit( wij) );
        ps3.push( ~Lit(I.eqi(s, i)));
        ps3.push( ~Lit(X[i-offset].eqi(s, j)));
        s.addClause(ps3);
      }
    }

  // -R=j, -w_1j ..., -wnj
  for(int j = rmin; j <= rmax; ++j) {
    if( !R.indomain(s, j) ) continue;
    ps.clear();
    ps.push( ~Lit(R.eqi(s, j)) );
    for(int i = imin; i <= imax; ++i) {
      Var wij = wfirst + idx(i, imin, j, rmin, n);
      ps.push( ~Lit(wij) );
    }
    s.addClause(ps);
  }
}

class cons_table;

/* Global constraints */

/* All different: each variable gets a distinct value */
class cons_alldiff : public cons
{
  static const int idx_undef;

  vector<cspvar> _x;
  bool _gac;
  // the matching, kept and updated between calls. No need to update
  // it on backtracking, whatever we had before is still valid
  vector<int> matching;    // from var to val
  vector<int> revmatching; // from val to var
  vector<int> matching0; // double buffering
  vector<int> revmatching0;
  vec<Lit> reason;

  int umin, umax; // the minimum and maximum values in the universe

  unsigned nmatched;

  typedef pair<bool, int> vertex;

  // all the following are here to avoid allocations
  vector<unsigned char> valvisited, varvisited; // for the BFS in matching,
                                       // whether it is in the stack in
                                       // tarjan's scc
  std::vector<int> varfrontier; // the frontier of the BFS
  std::vector<int> valfrontier;
  vector<int> valbackp, varbackp; // back pointers to reconstruct the
                                  // augmenting path
  vector<int> valvisited_toclear;
  vector<int> varvisited_toclear;

  vector<int> varindex, varlowlink;
  vector<int> valindex, vallowlink;

  vector<unsigned char> varhasfree; // did we reach a free value in the DFS?
  vector<unsigned char> valhasfree;

  vector<int> varcomp, valcomp;
  vector< vertex > components;
  vector< size_t > comp_limit;
  vector< bool > hallcomp; // is the scc a Hall set?
  vector< vertex > tarjan_stack;

  // SCC based decomposition
  vector<size_t> sccs;
  vector<size_t> scc_index;
  btptr scc_splitpoint_ptr; // a bool[]. Note that unlike the Gent et
                            // al paper, splitpoint[i] is true iff an
                            // scc starts at i

  vector<unsigned> scc_wake;
  vector<unsigned> scc_wake_idx;
  btptr scc_wake_size_ptr; // this is a single integer that is always
                           // 0 after propagate(). So if wake()
                           // records some changes and we backtrack
                           // before propagate() runs, scc_wake is
                           // automatically cleared

  vec<Lit>* reasons;

  bool varfree(size_t var) {
    return matching[var] < umin;
  }
  bool valfree(int val) {
    return revmatching[val-umin] < 0;
  }

  bool find_initial_matching(Solver& s) {
    nmatched = 0;
    greedy_matching(s);
    return find_matching(s);
  }

  // matching
  void greedy_matching(Solver& s);
  void clear_visited();
  bool find_matching(Solver& s);

  /* SCCs */
  // produce a clause that describes why there is no matching
  void explain_conflict(Solver& s, vec<Lit> &ps);
  // add to ps the reason that q is inconsistent. May require other
  // values to be explained and these are pushed to to_explain.
  void explain_value(Solver &s, int q,
                     vec<Lit>& ps,
                     vector<unsigned char>& explained,
                     vector<int>& to_explain);

  void tarjan_unroll_stack(Solver &s, vertex root);
  // start the dfs from variable var
  void tarjan_dfs(Solver& s, size_t var, size_t& index, bool conflict);
  // reset all structures that were touched by tarjan_*
  void tarjan_clear();
public:
  cons_alldiff(Solver &s, vector<cspvar> const& x, bool gac = true) :
    _x(x), _gac(gac)
  {
    set_priority(3);
    umin = _x[0].min(s);
    umax = _x[0].max(s);
    for(size_t i = 0; i != _x.size(); ++i) {
      umin = std::min( umin, _x[i].min(s) );
      umax = std::max( umax, _x[i].max(s) );
      s.schedule_on_dom(_x[i], this);
      s.wake_on_dom(_x[i], this, reinterpret_cast<void*>(i+1));
      s.wake_on_fix(_x[i], this);
    }

    matching.resize(_x.size(), umin-1);
    revmatching.resize(umax-umin+1, -1);

    varvisited.resize(_x.size(), false);
    varbackp.resize(_x.size(), umin-1);
    varindex.resize(_x.size(), idx_undef);
    varlowlink.resize(_x.size(), idx_undef);
    varcomp.resize(_x.size(), idx_undef );
    varhasfree.resize(_x.size(), false );

    valvisited.resize(umax-umin+1, false);
    valbackp.resize(umax-umin+1, -1);
    valindex.resize(umax-umin+1, idx_undef);
    vallowlink.resize(umax-umin+1, idx_undef);
    valcomp.resize(umax-umin+1, idx_undef);
    valhasfree.resize(umax-umin+1, false );

    comp_limit.push_back(0);

    reasons = new vec<Lit>[umax-umin+1 + _x.size()];

    sccs.resize(_x.size());
    scc_index.resize(_x.size());
    scc_splitpoint_ptr = s.alloc_backtrackable((_x.size()+1)*sizeof(bool));
    bool *scc_splitpoint = s.deref_array<bool>(scc_splitpoint_ptr);
    for(size_t i = 0; i != _x.size(); ++i) {
      sccs[i] = i;
      scc_index[i] = i;
      scc_splitpoint[i] = false;
    }
    scc_splitpoint[_x.size()] = true;

    scc_wake.resize(_x.size());
    scc_wake_idx.resize(_x.size());
    for(size_t i = 0; i != _x.size(); ++i) {
      scc_wake[i] = i;
      scc_wake_idx[i] = i;
    }
    scc_wake_size_ptr = s.alloc_backtrackable(sizeof(size_t));
    s.deref<size_t>(scc_wake_size_ptr) = 0;

    if( !find_initial_matching(s) ) {
      delete[] reasons;
      throw unsat();
    }
  }
  ~cons_alldiff() { delete[] reasons; }

  Clause* wake_advised(Solver &s, Lit p, void *);
  Clause* wake(Solver &s, Lit p);
  Clause *propagate(Solver& s);
  void clone(Solver& other);
  ostream& print(Solver &s, ostream& os) const;
  ostream& printstate(Solver& s, ostream& os) const;
};

const int cons_alldiff::idx_undef = std::numeric_limits<int>::max();

void cons_alldiff::clone(Solver& other)
{
  cons *con = new cons_alldiff(other, _x, _gac);
  other.addConstraint(con);
}

ostream& cons_alldiff::print(Solver& s, ostream& os) const
{
  os << "alldifferent(";
  for(size_t i = 0; i != _x.size(); ++i) {
    if( i ) os << ", ";
    os << cspvar_printer(s, _x[i]);
  }
  os << ")";
  return os;
}

ostream& cons_alldiff::printstate(Solver&s, ostream& os) const
{
  print(s, os);
  os << " (with ";
  for(size_t i = 0; i != _x.size(); ++i) {
    if( i ) os << ", ";
    os << cspvar_printer(s, _x[i]) << " in " << domain_as_set(s, _x[i]);
  }
  os << ")";
  return os;
}

void cons_alldiff::greedy_matching(Solver &s)
{
  const size_t n = _x.size();
  for(size_t i = 0; i != n; ++i) {
    for(int q = _x[i].min(s); q <= _x[i].max(s); ++q) {
      if( valfree(q) ) {
        matching[i] = q;
        revmatching[q-umin] = i;
        ++nmatched;
        break;
      }
    }
  }
}

void cons_alldiff::clear_visited()
{
  for(size_t i = 0; i != valvisited_toclear.size(); ++i)
    valvisited[ valvisited_toclear[i]-umin ] = false;
  for(size_t i = 0; i != varvisited_toclear.size(); ++i)
    varvisited[ varvisited_toclear[i] ] = false;
  valvisited_toclear.clear();
  varvisited_toclear.clear();
}

bool cons_alldiff::find_matching(Solver &s)
{
  const size_t n = _x.size();
  // find a free variable
  for(size_t fvar = 0; fvar != n && nmatched < n; ++fvar ) {
    if ( !varfree(fvar) ) continue;

    // do a bfs for an augmenting path
    varfrontier.push_back(fvar);
    int pathend = umin-1;
    do {
      while( !varfrontier.empty() ) {
        size_t var = varfrontier.back();
        varfrontier.pop_back();
        int qend = _x[var].max(s);
        for(int q = _x[var].min(s); q <= qend ; ++q) {
          if( !_x[var].indomainUnsafe(s, q) ) continue;
          if( matching[var] == q )
            continue;
          if( valvisited[q-umin] )
            continue;
          valbackp[q-umin] = var;
          if( valfree(q) ) {
            pathend = q;
            goto augment_path;
          }
          valfrontier.push_back(q);
          valvisited[q-umin] = true;
          valvisited_toclear.push_back(q);
        }
      }
      while( !valfrontier.empty() ) {
        int q = valfrontier.back();
        valfrontier.pop_back();
        int var = revmatching[q - umin];
        if( var < 0 || varvisited[var] ) continue;
        varfrontier.push_back(var);
        varbackp[var] = q;
        varvisited[var] = true;
        varvisited_toclear.push_back(var);
      }
    } while( !varfrontier.empty() );
    clear_visited();
    return false;
  augment_path:
    int val = pathend;
    int var = valbackp[pathend-umin];
    assert( var >= 0 );
    matching[var] = pathend;
    revmatching[pathend-umin] = var;
    while( (unsigned)var != fvar ) {
      val = varbackp[var];
      var = valbackp[val-umin];
      matching[var] = val;
      revmatching[val-umin] = var;
    }
    ++nmatched;
    valfrontier.clear();
    varfrontier.clear();
    clear_visited();
  }
  for(size_t i = 0; i != n; ++i)
    assert(!varfree(i));
  assert(nmatched == n);
  return true;
}

void cons_alldiff::explain_value(Solver &s, int q,
                                 vec<Lit>& ps,
                                 vector<unsigned char>& explained,
                                 vector<int>& to_explain)
{
  vector<unsigned char> hallvals( umax-umin+1, false );
  int scc = valcomp[q - umin];
  int sccsize = comp_limit[scc+1]-comp_limit[scc];
  if( sccsize == 1 ) {
    assert(!valhasfree[q-umin]);
    assert(revmatching[q-umin] >= 0);
    hallvals[q-umin] = true;
    explained[q-umin] = 2;
    scc = varcomp[ revmatching[q-umin] ];
  } else {
    // first gather the values of the Hall set
    for(size_t i = comp_limit[scc]; i != comp_limit[scc+1]; ++i) {
      vertex vx = components[i];
      if( vx.first ) continue; // var
      hallvals[vx.second-umin] = true;
      explained[vx.second-umin] = 2;
    }
  }
  // now describe that the variables in the SCC have been reduced to
  // form a Hall set
  for(size_t i = comp_limit[scc]; i != comp_limit[scc+1]; ++i) {
    vertex vx = components[i];
    if( !vx.first ) continue; // val
    cspvar v2 = _x[vx.second];
    pushifdef( ps, v2.r_min(s) );
    pushifdef( ps, v2.r_max(s) );
    for(int p = v2.min(s), pend = v2.max(s); p <= pend; ++p)
      if( !hallvals[p - umin] ) {
        if( explained[p-umin] )
          continue;
        if( v2.indomainUnsafe(s, p) ) {
          to_explain.push_back(p);
          explained[p-umin] = 1;
        }
        else
          ps.push( v2.r_neq(s, p) );
      }
  }
}

void cons_alldiff::explain_conflict(Solver &s, vec<Lit>& ps)
{
  // find a free var...
  size_t fvar;
  for(fvar = 0; !varfree(fvar); ++fvar)
    ;
  // ... and all SCCs reachable from it ...
  size_t index = 0;
  tarjan_dfs(s, fvar, index, true);
  cspvar v = _x[fvar];

  ps.clear();
  // describe the domain of fvar
  pushifdef( ps, v.r_min(s) );
  pushifdef( ps, v.r_max(s) );
  for(int q = v.min(s)+1; q < v.max(s); ++q)
    if( !v.indomainUnsafe(s, q) )
      ps.push( v.r_neq(s, q) );

  // and all the Hall sets that remove the rest of the values
  vector<unsigned char> explained( umax-umin+1, false );
  vector<int> to_explain;
  for(int q = v.min(s); q <= v.max(s); ++q) {
    if( !v.indomainUnsafe(s, q) ) continue;
    if( explained[q-umin] == 2 ) continue;
    explain_value(s, q, ps, explained, to_explain);
  }

  // note in this loop that to_explain.size() might change in the body
  // of the loop
  for(size_t i = 0; i != to_explain.size(); ++i) {
    if( explained[to_explain[i]-umin] == 2 )
      continue;
    explain_value(s, to_explain[i], ps, explained, to_explain);
  }

  tarjan_clear();
}

void cons_alldiff::tarjan_unroll_stack(Solver &s, vertex root)
{
  int scc = comp_limit.size()-1;
  int numvars = 0;
  size_t minvaridx = _x.size();
  hallcomp.push_back(true);
  vertex vp;
  do {
    vp = tarjan_stack.back(); tarjan_stack.pop_back();
    components.push_back(vp);
    if( vp.first ) {
      ++numvars;
      varcomp[vp.second] = scc;
      varvisited[vp.second] = false;
      if( varhasfree[vp.second] )
        hallcomp[scc] = false;
      minvaridx = std::min(minvaridx, scc_index[vp.second]);
    } else {
      valcomp[vp.second-umin] = scc;
      valvisited[vp.second-umin] = false;
      if( valhasfree[vp.second-umin] )
        hallcomp[scc] = false;
    }
  } while( vp != root );
  comp_limit.push_back(components.size());
  if( hallcomp[scc] ) {
    // split the constraint
    // find the original scc...
    if( numvars == 0 ) {
      int val = components[comp_limit[scc]].second;
      int var = revmatching[val-umin];
      assert( var >= 0 );
      minvaridx = scc_index[var];
    }
    unsigned vscc_start = minvaridx;
    bool *scc_splitpoint = s.deref_array<bool>(scc_splitpoint_ptr);
    for(; vscc_start > 0 && !scc_splitpoint[vscc_start]; --vscc_start)
      ;
    // ... now move this scc's variable to the beginning ...
    int moved = 0;
    unsigned i1 = vscc_start, i2 = vscc_start;
    for(; moved < numvars; ) {
      if( varcomp[ sccs[i1] ] == scc ) {
        if( i1 > i2 ) {
          std::swap(scc_index[sccs[i1]], scc_index[sccs[i2]]);
          std::swap(sccs[i1], sccs[i2]);
        }
        ++i1; ++i2;
        ++moved;
      } else
        ++i1;
    }
    // ... and mark the split.
    scc_splitpoint[i2] = true;
    return;
  }
  // not a hall set, so we can reach a free value
  for(size_t j = comp_limit[scc]; j != comp_limit[scc+1]; ++j) {
    vertex vx = components[j];
    if( vx.first )
      varhasfree[vx.second] = true;
    else
      valhasfree[vx.second-umin] = true;
  }
}

void cons_alldiff::tarjan_clear()
{
  for(size_t i = 0; i != varvisited_toclear.size(); ++i) {
    int var = varvisited_toclear[i];
    varindex[var] = idx_undef;
    varlowlink[var] = idx_undef;
    varcomp[var] = idx_undef;
    varhasfree[var] = false;
  }
  for(size_t i = 0; i != valvisited_toclear.size(); ++i) {
    int val = valvisited_toclear[i];
    valindex[val - umin] = idx_undef;
    vallowlink[val - umin] = idx_undef;
    valcomp[val - umin] = idx_undef;
    valhasfree[val - umin] = false;
  }
  for(size_t i = 0; i != comp_limit.size()-1; ++i)
    reasons[i].clear();
  valvisited_toclear.clear();
  varvisited_toclear.clear();
  components.clear();
  comp_limit.resize(1);
  hallcomp.clear();
}

void cons_alldiff::tarjan_dfs(Solver &s, size_t var, size_t& index,
                              bool conflict)
{
  varvisited_toclear.push_back(var);
  varindex[var] = index;
  varlowlink[var] = index;
  varvisited[var] = true;
  ++index;
  tarjan_stack.push_back( make_pair(true, var) );
  for(int q = _x[var].min(s), qend = _x[var].max(s); q <= qend; ++q) {
    if( !_x[var].indomainUnsafe(s, q) ) // no edge at all
      continue;
    if( matching[var] == q ) // edge is (q, var), not (var, q)
      continue;
    if( valindex[q-umin] == idx_undef  ) {
      int val = q;
      valvisited_toclear.push_back(val);
      valindex[val-umin] = index;
      vallowlink[val-umin] = index;
      valvisited[val-umin] = true;
      ++index;
      tarjan_stack.push_back( make_pair(false, val) );
      int var2 = revmatching[val-umin];
      if( var2 >= 0 ) {
        if( varindex[var2] == idx_undef  ) {
          tarjan_dfs(s, var2, index, conflict);
          vallowlink[val-umin] = std::min(vallowlink[val-umin], varlowlink[var2]);
        } else if( varvisited[var2] ) { // var2 is in stack
          vallowlink[val-umin] = std::min(vallowlink[val-umin], varindex[var2] );
        }
        valhasfree[val-umin] = varhasfree[var2];
      } else
        valhasfree[val-umin] = true;
      if( vallowlink[val-umin] == valindex[val-umin] ) {
        tarjan_unroll_stack( s, make_pair(false, val) );
      }
      varlowlink[var] = std::min(varlowlink[var], vallowlink[q-umin]);
    } else if( valvisited[q-umin] ) // q is in stack
      varlowlink[var] = std::min(varlowlink[var], valindex[q-umin] );
    if( valfree( q ) || valhasfree[q-umin] )
      varhasfree[var] = true;
    if( !conflict ) {
      int scc = valcomp[q-umin];
      if( scc == idx_undef ) continue; // q and var are in the same scc
      if( hallcomp[scc] ) { // Hall set
        if( reasons[scc].size() == 0 ) {
          vector<int> to_explain;
          vector<unsigned char> explained(umax-umin+1, false);
          explain_value(s, q, reasons[scc], explained, to_explain);
          assert(to_explain.empty());
        }
        _x[var].removef(s, q, reasons[scc]);
      }
    }
  }
  if( varlowlink[var] == varindex[var] ) {
    tarjan_unroll_stack( s, make_pair(true, var) );
  }
}

Clause* cons_alldiff::wake_advised(Solver &s, Lit p, void *advice)
{
  domevent de = s.event(p);

  assert( de.type == domevent::NEQ );

  size_t idx = reinterpret_cast<size_t>(advice)-1;
  size_t & scc_wake_size = s.deref<size_t>(scc_wake_size_ptr);
  if( scc_wake_idx[idx] >= scc_wake_size ) {
    unsigned prev = scc_wake[scc_wake_size];
    std::swap(scc_wake[scc_wake_size], scc_wake[scc_wake_idx[idx]]);
    std::swap(scc_wake_idx[idx], scc_wake_idx[prev]);
    ++scc_wake_size;
  }
  return 0L;
}

Clause* cons_alldiff::wake(Solver &s, Lit p)
{
  const size_t n = _x.size();
  domevent de = s.event(p);

  assert( de.type == domevent::EQ );

  reason.clear();
  reason.push(~p);
  for(size_t i = 0; i != n; ++i) {
    if( de.x == _x[i] ) continue;
    if( !_x[i].indomain(s, de.d) ) continue;
    DO_OR_RETURN(_x[i].removef(s, de.d, reason));
  }
  return 0L;
}

Clause* cons_alldiff::propagate(Solver &s)
{
  const size_t n = _x.size();
  assert(nmatched == n);
  bool valid = true; // is the current matching still valid?

  size_t & scc_wake_size = s.deref<size_t>(scc_wake_size_ptr);
  bool *scc_splitpoint = s.deref_array<bool>(scc_splitpoint_ptr);

  // the start index of those sccs that have been touched
  vector<unsigned> touch_ccs(n);
  vector<bool> touch_ccs_in(n);
  unsigned touch_ccs_size = 0;

  for(size_t i = 0; i != n; ++i) touch_ccs_in[i] = false;

  for(size_t w = 0; w != scc_wake_size; ++w) {
    size_t i = scc_wake[w];
    assert( !varfree(i) );
    if( !_x[i].indomainUnsafe( s, matching[i] ) ) {
      if(valid) {
        valid = false;
        revmatching0 = revmatching;
        matching0 = matching; // if we detect failure, we will copy
                              // matching0 back to matching so on
                              // backtracking we will still have a
                              // valid matching
      }
      revmatching[ matching[i]-umin ] = -1;
      matching[i] = umin-1;
      --nmatched;
    }
    if( !_gac ) continue;
    // we mark every variable we see so there is no danger of finding
    // the start of the scc potentially linear number of times. this
    // would be quadratic here.
    unsigned vscc_start = scc_index[i];
    for(; vscc_start > 0 && !scc_splitpoint[vscc_start]
          && !touch_ccs_in[vscc_start]; --vscc_start)
      touch_ccs_in[vscc_start] = true;
    if( !touch_ccs_in[vscc_start] ) {
      touch_ccs_in[vscc_start] = true;
      touch_ccs[touch_ccs_size] = vscc_start;
      ++touch_ccs_size;
    } // else we stopped on touch_ccs_in[..] being true, so we have
      // already added this scc to touch_ccs
  }
  scc_wake_size = 0;

  if( valid && !_gac ) {
    return 0L;
  }

  if( !valid && !find_matching(s) ) {
    reason.clear();
    explain_conflict(s, reason);
    Clause *r = Clause_new(reason);
    s.addInactiveClause(r);
    matching = matching0;
    revmatching = revmatching0;
    nmatched = n;
    return r;
  }

  assert(nmatched == n);

  if( !_gac )
    return 0L;

  size_t index = 0;
  for( unsigned scc = 0; scc != touch_ccs_size; ++scc ) {
    size_t idx = touch_ccs[scc];
    assert( varindex[ sccs[idx] ] == idx_undef );
    size_t eidx = idx;
    for(++eidx; !scc_splitpoint[eidx]; ++eidx)
      ;
    tarjan_dfs(s, sccs[idx], index, false);
    for(++idx; idx != eidx; ++idx) {
      if( varindex[ sccs[idx] ] == idx_undef )
        tarjan_dfs(s, sccs[idx], index, false);
    }
  }

  tarjan_clear();

  assert(nmatched == n);

  return 0L;
}

void post_alldiff(Solver &s, std::vector<cspvar> const &x, bool gac)
{
  vector<cspvar> xp;
  vector<int> rem;
  for(size_t i = 0; i != x.size(); ++i) {
    if( x[i].min(s) == x[i].max(s) ) {
      rem.push_back(x[i].min(s));
    } else
      xp.push_back(x[i]);
  }
  for(size_t i = 0; i != xp.size(); ++i) {
    for(size_t q = 0; q != rem.size(); ++q)
      xp[i].remove(s, q, NO_REASON);
  }
  if( xp.empty() ) return;

  cons *con = new cons_alldiff(s, xp, gac);
  s.addConstraint(con);
}

class cons_gcc;

/* Regular */
namespace regular {
  struct layered_fa {
    vector<transition> d;     // all transitions, ordered by layer
    vector<int> layer_trans;  // an index into the beginning of each
                              // layer in d
    vector<int> layer_states; // the first state of each layer
    vector<int> state_trans;  // index of first transition of each state
    vector<int> state_layer;  // the layer of each state
    int q0;
    set<int> F;

    int nlayers() const { return layer_states.size()-1; }
  };

  struct rejecting {
    typedef bool value_type;
    bool operator()(transition const& t) {
      return t.q1 == 0;
    }
  };

  struct compare_transition {
    typedef bool value_type;
    bool operator()(transition t1, transition t2) {
      // (t1.q0,t1.s,t1.q1) <=lex (t2.q0,t2.s, t2.q1)
      return t1.q0 < t2.q0 ||
        (t1.q0 == t2.q0 && (t1.s < t2.s ||
                            (t1.s == t2.s && t1.q1 <= t2.q1)));
    }
  };

  void unfold(automaton const& aut,
              Solver& solver,
              vector<cspvar> X,
              layered_fa& l)
  {
    // note for a string of length $n$ we need $n+1$ layers
    size_t n = X.size();
    using std::max;
    size_t ns = 0;
    for(size_t i = 0; i != aut.d.size(); ++i) {
      transition const& t = aut.d[i];
      ns = max(ns, max(t.q0, t.q1));
    }

    vector<transition> d = aut.d;
    sort(d.begin(), d.end(), compare_transition());

    // transitions and pointers into the transition table
    for(size_t i = 0; i != n; ++i) {
      // special handling for state 0. we assume it is at layer 0
      if( i == 0 )
        l.layer_states.push_back(0);
      else
        l.layer_states.push_back(ns*i+1);
      l.layer_trans.push_back(l.d.size());
      int prevq = -1;
      transition prevt(-1,-1,-1);
      for(size_t j = 0; j != aut.d.size(); ++j) {
        transition const& t = aut.d[j];
        if( t.q0 == 0 || t.q1 == 0 ) continue;
        if( !X[i].indomain(solver, t.s) ) continue;
        if( prevq < 0 || t.q0 != prevt.q0 ) {
          while( l.state_trans.size() < (unsigned)t.q0+i*ns ) {
            l.state_trans.push_back(l.d.size() );
            l.state_layer.push_back( i );
          }
          l.state_trans.push_back( l.d.size() );
          l.state_layer.push_back( i );
          prevq = t.q0;
        }
        l.d.push_back( transition(t.q0+i*ns, t.s, t.q1+(i+1)*ns) );
        prevt = t;
      }
    }
    // the last layer has all the states but no transitions
    l.layer_states.push_back(ns*n+1);
    l.layer_states.push_back(ns*(n+1)+1);
    l.layer_trans.push_back(l.d.size());
    l.layer_trans.push_back(l.d.size());

    while( l.state_trans.size() <= (n+1)*ns+1 ) {
      l.state_trans.push_back(l.d.size());
      l.state_layer.push_back(l.nlayers()-1);
    }

    // starting/accepting states
    l.q0 = aut.q0;
    for( set<int>::const_iterator i = aut.F.begin(), iend = aut.F.end();
         i != iend; ++i) {
      l.F.insert( *i + n*ns );
    }
  }

  void minimize_internal(layered_fa& l)
  {
    using std::max;
    size_t ns = 0;
    for(size_t i = 0; i != l.d.size(); ++i) {
      transition& t = l.d[i];
      ns = max(ns, max(t.q0, t.q1));
    }

    size_t accepting=0; // all accepting states are merged into one
    vector<int> remap(ns+1);
    // map from set of (s, q1) to all the q0s that have this. In the
    // end we merge all q0s in the same bucket. inefficient, this
    // should be a trie
    typedef std::map< set< pair<int, size_t> >, vector< size_t > > mergemap;
    mergemap tomerge;
    for(int n = l.nlayers()-1; n >= 0; --n) {
      tomerge.clear();
      for(int s = l.layer_states[n]; s != l.layer_states[n+1]; ++s) {
        remap[s] = s;
        set< pair<int, size_t> > out;
        for(int i = l.state_trans[s]; i != l.state_trans[s+1]; ++i) {
          // note we both gather the outgoing (s,q1) of each q0 AND
          // update the transitions to reflect the states we have
          // fixed already in layer n+1
          transition& t = l.d[i];
          t.q1 = remap[t.q1];
          if( t.q1 != 0 )
            out.insert( make_pair( t.s, t.q1 ) );
        }
        if( out.empty() ) { // final layer state or no path to accepting
          if( l.F.find(s) == l.F.end() ) { // final layer rejecting or
                                           // non-final layer and no
                                           // path to accepting
            remap[s] = 0;
            continue;
          } // else final layer accepting state
        }
        tomerge[out].push_back(s);
      }
      for(mergemap::const_iterator i=tomerge.begin(), iend=tomerge.end();
          i != iend; ++i) {
        vector<size_t> const & mergestates = i->second;
        for(vector<size_t>::const_iterator j = mergestates.begin(),
              jend = mergestates.end();
            j != jend; ++j) {
          if( i->first.empty() ) // last layer, therefore accepting
            accepting = mergestates[0];
          remap[*j] = mergestates[0];
        }
      }
    }
    assert( accepting > 0 );
    l.F.clear();
    l.F.insert( accepting );
  }

  void minimize_unfolded(layered_fa& l)
  {
    // FIXME: for now this only minimizes DFAs, but if the automaton
    // is non-deterministic, we should do heuristic minimization as in
    // [KNW CPAIOR09]
    minimize_internal(l);
  }

  void gather_reachable(layered_fa& aut, set<size_t>& r)
  {
    using std::max;
    size_t ns = 0;
    for(size_t i = 0; i != aut.d.size(); ++i) {
      transition& t = aut.d[i];
      ns = max(ns, max(t.q0, t.q1));
    }
    vector<int> remap(ns, 0);

    std::list<size_t> Q;
    Q.push_back(aut.q0);
    r.insert(aut.q0);
    while( !Q.empty() ) {
      size_t s = Q.front();
      Q.pop_front();

      for(int t = aut.state_trans[s]; t != aut.state_trans[s+1]; ++t) {
        transition const& tr = aut.d[t];
        assert(tr.q0 == s);
        if( tr.q1 == 0 ) continue;
        if( r.find(tr.q1) == r.end() ) {
          r.insert(tr.q1);
          Q.push_back(tr.q1);
        }
      }
    }
  }

  void ensure_var(Solver& s, Var& var)
  {
    if( var == var_Undef )
      var = s.newVar();
  }
};

void post_regular(Solver &s, vector< cspvar > const& vars,
                  regular::automaton const& aut,
                  bool gac)
{
  using namespace regular;
  layered_fa l;
  unfold(aut, s, vars, l);
  minimize_unfolded(l);

  assert(l.F.size() == 1);

  set<size_t> r;
  gather_reachable(l, r);

  if( r.find(*l.F.begin()) == r.end() )
    throw unsat();

  vector<Var> sv(l.state_trans.size(), var_Undef); // State variables
  vector<Var> tv(l.d.size(), var_Undef); // Transition variables

  vec<Lit> ps, ps2;

  ensure_var(s, sv[l.q0]);
  ps.push( Lit(sv[l.q0]) ); // Must start in the initial state
  s.addClause(ps);

  ps.clear();
  ensure_var(s, sv[*l.F.begin()]);
  ps.push( Lit(sv[*l.F.begin()]) ); // Must end in the accepting state
  s.addClause(ps);

  typedef map< int, vector<int> > incoming_map_t;
  incoming_map_t imap; // Map each state to a vector of states
                       // with outgoing edges to it
  typedef map< pair<int, int>, vector<int> > trans_map_t;
  trans_map_t tmap; // Map each V=d to the vector of transitions

  for( set<size_t>::iterator si = r.begin(), siend = r.end();
       si != siend; ++si) {
    size_t st = *si;
    size_t layer = l.state_layer[st];
    ensure_var(s, sv[st]);

    if( (unsigned)l.state_layer[st] == vars.size() ) continue;

    ps2.clear();
    ps2.push( ~Lit( sv[st] ) );
    for(int i = l.state_trans[st]; i != l.state_trans[st+1]; ++i) {
      transition tr = l.d[i]; // All outgoing transitions from `st'

      if( tr.q1 == 0 ) continue; // continue if there is no path to an accepting state
      assert( tr.q0 == st );
      ensure_var(s, tv[i]);
      ensure_var(s, sv[tr.q1]);

      imap[tr.q1].push_back(tv[i]);
      tmap[make_pair(layer, tr.s)].push_back(tv[i]);

      ps.clear();
      ps.push( ~Lit( tv[i] ) );
      ps.push( Lit( sv[tr.q0] ) );
      s.addClause(ps); // Transition from q0 to q1 implies that state is in the
                       // next state (3b in [Bacchus, CP07])

      ps.clear();
      ps.push( ~Lit( tv[i] ) );
      ps.push( Lit( sv[tr.q1] ) );
      s.addClause(ps); // Transition from q0 to q1 implies that state was in the
                       // previous step (3a in [Bacchus, CP07])
      ps.clear();
      ps.push( ~Lit( tv[i] ) );
      ps.push( Lit( vars[layer].eqi(s, tr.s) ) ); // Transition happened means
                                                  // corresponding variable took that
                                                  // value (3c in [Bacchus, CP07])
      s.addClause(ps);

      ps2.push( Lit( tv[i] ) );
    }
    s.addClause(ps2); // At least one outgoing transition is true
                      // (4a in [Bacchus, CP07])
    ps2.clear();
  }

  if( !gac ) return;

  // For each state, at least one incoming transition must be true
  // (4b in [Bacchus, CP07])
  for (incoming_map_t::const_iterator mi = imap.begin(), miend = imap.end();
       mi != miend; ++mi) {
    int q1 = mi->first;
    const vector<int> &incoming = mi->second;

    ps2.push(~Lit(sv[q1]));
    for (vector<int>::const_iterator vi = incoming.begin(), viend = incoming.end();
         vi != viend; ++vi) {
      Var tv = *vi;
      ps2.push(Lit(tv));
    }

    s.addClause(ps2);
    ps2.clear();
  }

  // For each V=d, at least one supporting transition must be enabled
  // (5 in [Bacchus, CP07])
  for(int i = 0; i != int(vars.size()); ++i) {
      cspvar var = vars[i];
      for(int val = var.min(s); val <= var.max(s); ++val) {
          if( !var.indomain(s, val) ) continue;
          const vector<int> &transitions = tmap[make_pair(i, val)];

          ps2.push(~Lit(var.eqi(s, val)));

          for (vector<int>::const_iterator ti = transitions.begin(), tiend = transitions.end();
               ti != tiend; ++ti) {
              Var tv = *ti;
              ps2.push(Lit(tv));
          }

          s.addClause(ps2);
          ps2.clear();
      }
  }
}

namespace cumulative {
  bool fixed(Solver& s, cspvar x)
  {
    return x.min(s) == x.max(s);
  }

  int value(Solver& s, cspvar x)
  {
    assert(fixed(s, x));
    return x.min(s);
  }

  typedef std::map<cspvar, vector<Lit> > running_map;
  std::map<Solver *, running_map> solver_running;
}

void post_cumulative(Solver& s, vector<cspvar> const& start,
                     vector<cspvar> const& dur,
                     vector<cspvar> const& req,
                     cspvar cap)
{
  using namespace cumulative;

  const size_t n = start.size();
  assert( dur.size() == n && req.size() == n );
  int mint = start[0].min(s);
  int maxt = start[0].max(s)+dur[0].max(s);
  for(size_t i = 1; i != n; ++i) {
    mint = std::min(mint, start[i].min(s));
    maxt = std::max(maxt, start[i].max(s) + dur[i].max(s));
  }

  Var truevar = s.newVar();
  Lit truelit = Lit(truevar);
  Lit falselit = ~Lit(truevar);
  s.uncheckedEnqueue(truelit);

  for(int t = mint; t != maxt+1; ++t) {
    vector<Lit> running(n, lit_Undef);
    bool unfixed_req = false;
    for(size_t i = 0; i != n; ++i) {
      if( start[i].min(s) > t || start[i].max(s) < t )
        continue;

      std::vector<Lit> & task_running = cumulative::solver_running[&s][start[i]];
      if( !task_running.empty() &&
          task_running[t - start[i].omin(s)] != lit_Undef) {
        running[i] = task_running[t - start[i].omin(s)];
        if( !fixed(s, req[i]) ) unfixed_req = true;
        continue;
      }
      if( task_running.empty() )
        task_running.resize( maxt-start[i].omin(s), lit_Undef );

      Lit started, finished;

      started = start[i].e_leq(s, t);
      if( fixed(s, dur[i]) ) {
        if( value(s, dur[i]) == 1 ) {
          running[i] = start[i].e_eq(s, t);
          task_running[t - start[i].omin(s)] = running[i];
          continue;
        }
        if( start[i].min(s)+value(s, dur[i]) > t )
          finished = falselit;
        else if( start[i].max(s)+value(s, dur[i]) <= t )
          finished = truelit;
        else
          finished = start[i].e_leq(s, t-value(s, dur[i]));
      } else {
        // s[i]+d[i] <= t <=> finished
        // excellent use for neg_re, if we had it
        cspvar b = s.newCSPVar(0, 1);
        vector<cspvar> v(2);
        vector<int> w(2);
        v[0] = start[i]; w[0] = 1;
        v[1] = dur[i];   w[1] = 1;
        post_lin_leq_iff_re(s, v, w, -t, b);
        finished = b.e_eq(s, 1);
      }


      running[i] = Lit(s.newVar());
      task_running[t - start[i].omin(s)] = running[i];
      vec<Lit> ps, ps2;
      ps.push( started );
      ps.push( ~running[i] );
      s.addClause(ps);
      ps.clear();

      ps.push( ~finished );
      ps.push( ~running[i] );
      s.addClause(ps);
      ps.clear();

      ps2.push(~started);
      ps2.push(finished);
      ps2.push( running[i] );
      s.addClause(ps2);

      if( !fixed(s, req[i]) ) unfixed_req = true;
    }

    if( unfixed_req ) { // at least one task has unfixed resource
                        // requirement
      assert(0);
    } else {
      vector<int> w;
      vector<int> v;
      for(size_t i = 0; i != n; ++i) {
        if( running[i] == lit_Undef ) continue;
        w.push_back( value(s, req[i]) );
        v.push_back( var(running[i]) );
      }
      cspvar tcap = s.newCSPVar(0, cap.max(s));
      post_leq(s, tcap, cap, 0);
      post_pb(s, v, w, 0, tcap);
    }
  }
}

namespace table {
  typedef map<int, int> state_trans;
  typedef map< int, state_trans > transmap;

  set<int> gather_sigma(Solver& s, vector<cspvar> const& x)
  {
    set<int> sigma;
    for(size_t i = 0; i != x.size(); ++i) {
      cspvar xi = x[i];
      for(int j = xi.min(s); j <= xi.max(s); ++j)
        if( xi.indomain(s, j) )
          sigma.insert(j);
    }
    return sigma;
  }

  void add_det_word(vector<int> const& tuple, transmap& d,
                    int& numstates,
                    const int q0, const int accepting)
  {
    size_t i = 0;
    int q = q0;
    while(i < tuple.size()-1) {
      state_trans& qtrans = d[q];
      state_trans::iterator sti = qtrans.find(tuple[i]);
      if(sti != qtrans.end() )
        q = sti->second;
      else {
        q = numstates;
        ++numstates;
        qtrans[tuple[i]] = q;
      }
      ++i;
    }
    d[q][tuple[i]] = accepting;
  }

  void add_absorbing_state(transmap& d, set<int> const& sigma,
                           const int absorbing)
  {
    for(transmap::iterator i = d.begin(), iend = d.end();
        i != iend; ++i) {
      state_trans& qtrans = i->second;
      for(set<int>::const_iterator j = sigma.begin(), jend = sigma.end();
          j != jend; ++j) {
        state_trans::iterator sti = qtrans.find(*j);
        if( sti == qtrans.end() )
          qtrans[*j] = absorbing;
      }
    }
  }

  void verify_automaton(vector< vector<int> > const& tuples,
                        transmap& a,
                        const int q0,
                        set<int> F)
  {
    for(size_t i = 0; i != tuples.size(); ++i) {
      vector<int> const & t = tuples[i];
      int q = q0;
      for(size_t j = 0; j != t.size(); ++j) {
        assert( a[q].find(t[j]) != a[q].end() );
        q = a[q][t[j]];
      }
      assert(F.find(q) != F.end());
    }
  }

  regular::automaton table_to_dfa(vector< vector<int> > const& tuples,
                                  set<int> const& sigma,
                                  bool positive)
  {
    int q0 = 1;
    set<int> aF;
    vector<regular::transition> ad;

    const int absorbing = 2;
    const int accepting = 3;
    int numstates = 4;

    transmap d;
    for(size_t i = 0; i != tuples.size(); ++i)
      add_det_word(tuples[i], d, numstates, q0, accepting);
    add_absorbing_state(d, sigma, absorbing);

    if( positive )
      aF.insert(accepting);
    else
      aF.insert(absorbing);

    verify_automaton(tuples, d, q0, aF);

    for( transmap::const_iterator i = d.begin(), iend = d.end();
         i != iend; ++i) {
      int q0 = i->first;
      for( state_trans::const_iterator
             j = i->second.begin(),
             jend = i->second.end();
           j != jend; ++j)
        ad.push_back(regular::transition( q0, j->first, j->second));
    }

    regular::automaton a(ad, q0, aF);
    return a;
  }

  // Fahiem's one-variable-per-tuple encoding. I don't know why I
  // called it ac4. Also enforces GAC for short tables (tables with
  // STAR_CONSTANT instead of a value in some tuples/positions)
  void post_positive_table_ac4(Solver &s, std::vector<cspvar> const &x,
                               std::vector<std::vector<int>> const &tuples) {
    vector< map<int, set<Var> > > valtuples(x.size());
    vec<Lit> ps;
    for(size_t i = 0; i != tuples.size(); ++i) {
      if(tuples[i].size() != x.size())
        throw non_table();
      Var is = s.newVar();
      for (size_t j = 0; j != tuples[i].size(); ++j) {
        if (tuples[i][j] == STAR_CONSTANT) {
          for (int q = x[j].min(s); q <= x[j].max(s); ++q)
            if (x[j].indomain(s, q))
              valtuples[j][q].insert(is);
          continue;
        }
        ps.clear();
        ps.push(~Lit(is));
        ps.push(x[j].e_eq(s, tuples[i][j]));
        s.addClause(ps);
        valtuples[j][tuples[i][j]].insert(is);
      }
    }

    for(size_t i = 0; i != x.size(); ++i) {
      for(int j = x[i].min(s); j <= x[i].max(s); ++j) {
        vec<Lit> ps;
        ps.push( x[i].e_neq( s, j ) );
        for(set<Var>::const_iterator si = valtuples[i][j].begin(),
              siend = valtuples[i][j].end(); si != siend; ++si)
          ps.push( Lit(*si) );
        s.addClause(ps);
      }
    }
  }

  void post_positive_table_regular(Solver &s, std::vector<cspvar> const& x,
                                   std::vector< std::vector<int> > const& tuples)
  {
    set<int> sigma = gather_sigma(s, x);
    regular::automaton a = table_to_dfa(tuples, sigma, true);
    post_regular(s, x, a);
  }
}


void post_positive_table(Solver &s, std::vector<cspvar> const& x,
                         std::vector< std::vector<int> > const& tuples)
{
  table::post_positive_table_ac4(s, x, tuples);
}

// Here we just post everything as a clause, better encodings later
void post_negative_table(Solver &s, std::vector<cspvar> const& x,
                         std::vector< std::vector<int> > const& tuples)
{
  for(size_t i = 0; i != tuples.size(); ++i) {
    if(tuples[i].size() != x.size())
      throw non_table();
    vec<Lit> ps;
    for(size_t j = 0; j != tuples[i].size(); ++j)
      if (tuples[i][j] != STAR_CONSTANT)
        ps.push( x[j].r_eq( s, tuples[i][j] ) );
    s.addClause(ps);
  }
}

void post_lex_common(Solver &s, std::vector<cspvar> const& x,
                     std::vector<cspvar> const& y,
                     std::vector<Var> const& b,
                     std::vector<Var> const& c)
{
    /* we post the decomposition
     * bi -> xi <= yi
     * ci <-> xi < yi
     * bi /\ -ci -> b(i+1)
     * b0
     *
     * Note that the first constraint is a one-way implication. If it
     * is equivalence we would need a third Boolean var to get the
     * same effect.
     */

    using std::size_t;
    size_t n = x.size();
    vec<Lit> ps;
    for(size_t i = 0; i != n; ++i) {
        post_leq_re_li(s, x[i], y[i], 0, Lit(b[i]));
        post_leq_re(s, x[i], y[i], -1, Lit(c[i]));

        if( i < n-1 ) {
            ps.clear();
            ps.push(~Lit(b[i]));
            ps.push(Lit(c[i]));
            ps.push(Lit(b[i+1]));
            s.addClause(ps);
        }
    }
    s.enqueue(Lit(b[0]));
}

void post_lex_leq(Solver &s, std::vector<cspvar> const& x,
                  std::vector<cspvar> const& y)
{
    assert(x.size() == y.size());
    std::size_t n = x.size();
    std::vector<Var> b(n), c(n);
    for(size_t i = 0; i != n; ++i) {
        b[i] = s.newVar();
        c[i] = s.newVar();
    }
    post_lex_common(s, x, y, b, c);
}

void post_lex_less(Solver &s, std::vector<cspvar> const& x,
                   std::vector<cspvar> const& y)
{
    assert(x.size() == y.size());
    size_t n = x.size();
    std::vector<Var> b(n), c(n);
    for(size_t i = 0; i != n; ++i) {
        b[i] = s.newVar();
        c[i] = s.newVar();
    }
    post_lex_common(s, x, y, b, c);
    vec<Lit> ps;
    ps.push(~Lit(b[n-1]));
    ps.push(Lit(c[n-1]));
    s.addClause(ps);
}

} // namespace minicsp
