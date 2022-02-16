/**************************************************************************
minicsp

Copyright 2009--2011 George Katsirelos

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

Parts of this file were distributed under the following license as
part of MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***************************************************************************/


#ifndef __MINICSP_SOLVERTYPES_HPP
#define __MINICSP_SOLVERTYPES_HPP

#include <cassert>
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <memory>

#include "minicsp/mtl/Alg.h"
#include "minicsp/mtl/Vec.h"

namespace minicsp {

class Solver;

// copied from boost/preprocessor/cat.hpp
#define BOOST_PP_CAT(a, b) BOOST_PP_CAT_I(a, b)
#define BOOST_PP_CAT_I(a, b) a ## b

// the number of queues. 0--3 because it can be encoded in the unused low
// order bits of a pointer in a system with word-alignment
#define MAX_PRIORITY 3

//=================================================================================================
// Variables, literals, lifted booleans, clauses:


// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
inline const Var var_Undef = -1;

class Lit {
    int     x;
 public:
    Lit() : x(2*var_Undef)                                              { }   // (lit_Undef)
    explicit Lit(Var var, bool sign = false) : x((var+var) + (int)sign) { }

    // Don't use these for constructing/deconstructing literals. Use the normal constructors instead.
    friend int  toInt       (Lit p);  // Guarantees small, positive integers suitable for array indexing.
    friend Lit  toLit       (int i);  // Inverse of 'toInt()'
    friend Lit  operator   ~(Lit p);
    friend bool sign        (Lit p);
    friend int  var         (Lit p);
    friend Lit  unsign      (Lit p);
    friend Lit  id          (Lit p, bool sgn);

    bool operator == (Lit p) const { return x == p.x; }
    bool operator != (Lit p) const { return x != p.x; }
    bool operator <  (Lit p) const { return x < p.x;  } // '<' guarantees that p, ~p are adjacent in the ordering.
};

inline  int  toInt       (Lit p)           { return p.x; }
inline  Lit  toLit       (int i)           { Lit p; p.x = i; return p; }
inline  Lit  operator   ~(Lit p)           { Lit q; q.x = p.x ^ 1; return q; }
inline  bool sign        (Lit p)           { return p.x & 1; }
inline  int  var         (Lit p)           { return p.x >> 1; }
inline  Lit  unsign      (Lit p)           { Lit q; q.x = p.x & ~1; return q; }
inline  Lit  id          (Lit p, bool sgn) { Lit q; q.x = p.x ^ (int)sgn; return q; }

const Lit lit_Undef(var_Undef, false);  // }- Useful special constants.
const Lit lit_Error(var_Undef, true );  // }


//=================================================================================================
// Lifted booleans:


class lbool {
    char     value;
    explicit lbool(int v) : value(v) { }

public:
    lbool()       : value(0) { }
    lbool(bool x) : value((int)x*2-1) { }
    int toInt(void) const { return value; }

    bool  operator == (lbool b) const { return value == b.value; }
    bool  operator != (lbool b) const { return value != b.value; }
    lbool operator ^ (bool b) const { return b ? lbool(-value) : lbool(value); }

    friend int   toInt  (lbool l);
    friend lbool toLbool(int   v);
};
inline int   toInt  (lbool l) { return l.toInt(); }
inline lbool toLbool(int   v) { return lbool(v);  }

const lbool l_True  = toLbool( 1);
const lbool l_False = toLbool(-1);
const lbool l_Undef = toLbool( 0);

// ======================================================================
// Pointer variant
//
// Store one of a set of pointers, using the free bits at the bottom
// to store the discriminator

namespace tagged_pointer_detail {

template <typename T> inline constexpr int tagged_free_bits()
{
    int x = alignof(T*);
    int bit{0};
    while(x>>=1)
      ++bit;
    return bit;
}

// tagged pointer: store a pointer to T and a tag of N bits. The tag
// can be accessed either as an integer (::tag()), or as a set of bits
// (::[g,s,res]et_bit()). It requires that the alignment requirements
// for T are strict enough to leave N bits free at the bottom
template <typename T, int N> struct tagged_pointer_N {
  static_assert(N <= tagged_free_bits<T>(),
                "pointer type for T does not have N free bits");
  // mask: all pointers satisfy p & ~mask == 0, if N is set to honor
  // the alignment of T
  static constexpr uintptr_t mask = ~uintptr_t(0) << N;
  uintptr_t p{0};

  void set(T *pointer, std::size_t tag = 0) {
    assert((tag & mask) == 0);
    auto pdata = reinterpret_cast<std::size_t>(pointer);
    assert((pdata & ~mask) == 0);
    p = pdata | tag;
  }
  T *get() { return reinterpret_cast<T *>(p & mask); }
  const std::size_t tag() { return p & ~mask; }
  template <int Bit> bool get_bit() {
    static_assert(Bit <= N, "This type only holds N bits");
    return p & ~(1u << (Bit - 1));
  }
  template <int Bit> void set_bit() {
    static_assert(Bit <= N, "This type only holds N bits");
    p |= (1u << (Bit - 1));
  }
  template <int Bit> void reset_bit() {
    static_assert(Bit <= N,  "This type only holds N bits");
    p &= ~(1u << (Bit - 1));
  }

};

template <typename T>
using tagged_pointer = tagged_pointer_N<T, tagged_free_bits<T>()>;

template <typename... Types> inline constexpr int tagged_free_bits_variant() {
  return std::min({tagged_free_bits<Types>()...});
}

template <typename T, typename... Types> inline constexpr int contained() {
  std::initializer_list<bool> b{false, std::is_same<T, Types>::value...};
  return std::max(b);
}

template <typename T, typename... Types> inline constexpr int indexof() {
  bool b[sizeof...(Types)]{std::is_same<T, Types>::value...};
  for (std::size_t i = 0; i != sizeof...(Types); ++i)
    if (b[i] == true)
      return i;
  return sizeof...(Types) + 1;
}

template <typename... Types> inline constexpr bool has_duplicates_inner() {
  return false;
}
template <typename T1, typename... Types>
inline constexpr bool has_duplicates_inner(T1 *, Types *... targs) {
  return contained<T1, Types...>() || has_duplicates_inner(targs...);
}

template <typename T1, typename... Types>
inline constexpr bool has_duplicates() {
  // we have to create actual (dummy) arguments, because we cannot
  // overload functions on the template arguments
  return has_duplicates_inner(static_cast<T1 *>(nullptr),
                              static_cast<Types *>(nullptr)...);
}

} //namespace tagged_pointer_detail

template <typename... Types> struct pointer_variant {
  static constexpr int freebits =
      tagged_pointer_detail::tagged_free_bits_variant<Types...>();
  static constexpr int maxtypes = 1 << (freebits - 1);
  static_assert(maxtypes > sizeof...(Types),
                "Not enough free bits to store this many types");
  static_assert(!tagged_pointer_detail::has_duplicates<Types...>(),
                "Cannot differentiate between pointers of the same type");

  // default constructor etc
  pointer_variant() = default;
  pointer_variant(pointer_variant<Types...> const &o) = default;
  template <typename T> pointer_variant(T *p) { set(p); }

  template <typename T> void set(T *p) {
    static_assert(tagged_pointer_detail::contained<T, Types...>(),
                  "T is not one of the pointer types stored");
    data.set(p, 1+ tagged_pointer_detail::indexof<T, Types...>());
  }
  void reset() { data.set(nullptr, 0); }
  template <typename T> T *get() {
    static_assert(tagged_pointer_detail::contained<T, Types...>(),
                  "T is not one of the pointer types stored");
    if (has<T>())
      return static_cast<T *>(data.get());
    else
      return nullptr;
  }

  operator bool() { return !null(); }

  // note we use 1 + index, because 0 is reserved for no type. Reduces
  // the number of usable types by 1, but allows has() to return the
  // right thing even for default constructed instances and to
  // differentiate between (T1*)nullptr and (T2*)nullptr
  template <typename T> bool has() {
    static_assert(tagged_pointer_detail::contained<T, Types...>(),
                  "T is not one of the pointer types stored");
    return (data.tag() == 1 + tagged_pointer_detail::indexof<T, Types...>());
  }

  bool empty() { return data.tag() == 0; }
  bool null() { return empty() || data.get() == nullptr; }

  tagged_pointer_detail::tagged_pointer_N<void, freebits> data;
};

//=================================================================================================
// Clause -- a simple class for representing a clause:

class Clause;

template<class V>
Clause* Clause_new(V&& ps, bool learnt = false,
                   Lit effect = lit_Undef);

class Clause {
    uint32_t size_etc;
    union { float act; uint32_t abst; } extra;
    Lit     data[0];

public:
    void calcAbstraction() {
        uint32_t abstraction = 0;
        for (int i = 0; i < size(); i++)
            abstraction |= 1 << (var(data[i]) & 31);
        extra.abst = abstraction;  }

    // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
    template<class V>
    Clause(V&& ps, bool learnt, Lit effect) {
        size_etc = (ps.size() << 3) | (uint32_t)learnt;
        for (int i = 0; i < ps.size(); i++) {
          assert(var(ps[i]) != var_Undef);
          data[i] = ps[i];
          if( data[i] == effect ) std::swap(data[0], data[i]);
        }
        if (learnt) extra.act = 0; else calcAbstraction(); }

    // -- use this function instead:
    template<class V>
    friend Clause* Clause_new(V&& ps, bool learnt, Lit effect) {
        assert(sizeof(Lit)      == sizeof(uint32_t));
        assert(sizeof(float)    == sizeof(uint32_t));
        void* mem = malloc(sizeof(Clause) + sizeof(uint32_t)*(ps.size()));
        return new (mem) Clause(ps, learnt, effect); }

    int          size        ()      const   { return size_etc >> 3; }
    void         shrink      (int i)         { assert(i <= size()); size_etc = (((size_etc >> 3) - i) << 3) | (size_etc & 7); }
    void         pop         ()              { shrink(1); }
    bool         learnt      ()      const   { return size_etc & 1; }
    uint32_t     mark        ()      const   { return (size_etc >> 1) & 3; }
    void         mark        (uint32_t m)    { size_etc = (size_etc & ~6) | ((m & 3) << 1); }
    const Lit&   last        ()      const   { return data[size()-1]; }

    // NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
    //       subsumption operations to behave correctly.
    Lit&         operator [] (int i)         { return data[i]; }
    Lit          operator [] (int i) const   { return data[i]; }
    operator const Lit* (void) const         { return data; }

    float&       activity    ()              { return extra.act; }
    uint32_t     abstraction () const { return extra.abst; }

    Lit          subsumes    (const Clause& other) const;
    void         strengthen  (Lit p);

    Lit const *begin() const { return static_cast<Lit const*>(*this); }
    Lit const *end() const { return begin() + size(); }
};

/*_________________________________________________________________________________________________
|
|  subsumes : (other : const Clause&)  ->  Lit
|
|  Description:
|       Checks if clause subsumes 'other', and at the same time, if it can be used to simplify 'other'
|       by subsumption resolution.
|
|    Result:
|       lit_Error  - No subsumption or simplification
|       lit_Undef  - Clause subsumes 'other'
|       p          - The literal p can be deleted from 'other'
|________________________________________________________________________________________________@*/
inline Lit Clause::subsumes(const Clause& other) const
{
    if (other.size() < size() || (extra.abst & ~other.extra.abst) != 0)
        return lit_Error;

    Lit        ret = lit_Undef;
    const Lit* c  = (const Lit*)(*this);
    const Lit* d  = (const Lit*)other;

    for (int i = 0; i < size(); i++) {
        // search for c[i] or ~c[i]
        for (int j = 0; j < other.size(); j++)
            if (c[i] == d[j])
                goto ok;
            else if (ret == lit_Undef && c[i] == ~d[j]){
                ret = c[i];
                goto ok;
            }

        // did not find it
        return lit_Error;
    ok:;
    }

    return ret;
}


inline void Clause::strengthen(Lit p)
{
    remove(*this, p);
    calcAbstraction();
}

/**********************
 *
 * CSP stuff
 *
 **********************/

class Solver;
class cons;

// a container for operations, lightweight enough to pass by value
class cspvar
{
  friend class Solver;
  int _id;
 public:
  cspvar() : _id(-1) {}
  explicit cspvar(int id) : _id(id) {}

  int id() const { return _id; }
  bool valid() const { return _id >= 0; }

  // current domain
  bool indomain(Solver& s, int d) const;
  // this performs no check whether ind(d)
  bool indomainUnsafe(Solver& s, int d) const;
  int min(Solver& s) const;
  int max(Solver& s) const;

  int domsize(Solver& s) const;

  // the original domain
  int omin(Solver& s) const;
  int omax(Solver& s) const;

  /* the reason can be an existing clause, a clause to be generated
     from a vec<Lit> or a cons, which will lazily be called on to
     generate a clause later, if needed.

     These return a Clause *, which is NULL if the operation does not
     cause a conflict, otherwise return the reason, which is a failed
     clause. Note if the reason is NULL (NO_REASON) and there is a
     failure, they throw unsat().
  */
  Clause *remove(Solver& s, int d, Clause *c);
  template<typename V>
  Clause *remove(Solver& s, int d, V& reason);
  Clause *remove(Solver& s, int d, cons *c);

  Clause *removerange(Solver &s, int rb, int re, Clause *c);
  template<typename V>
  Clause *removerange(Solver &s, int rb, int re, V& reason);
  Clause *removerange(Solver &s, int rb, int re, cons *c);

  Clause *setmin(Solver& s, int m, Clause *c);
  template<typename V>
  Clause *setmin(Solver& s, int m, V& reason);
  Clause *setmin(Solver& s, int m, cons *c);

  Clause *setmax(Solver& s, int m, Clause *c);
  template<typename V>
  Clause *setmax(Solver& s, int m, V& reason);
  Clause *setmax(Solver& s, int m, cons *c);

  Clause *assign(Solver& s, int d, Clause *c);
  template<typename V>
  Clause *assign(Solver& s, int d, V& reason);
  Clause *assign(Solver& s, int d, cons *c);

  /* the *f versions of the template<V> overloads above fill in the
     reason with the appropriate e_ literal. The literal is push()ed
     into reason and pop()ed back before returning, so reason is
     not const but is not modified.
  */
  template<typename V>
  Clause *removef(Solver& s, int d, V& reason);
  template<typename V>
  Clause *setminf(Solver& s, int m, V& reason);
  template<typename V>
  Clause *setmaxf(Solver& s, int m, V& reason);
  template<typename V>
  Clause *assignf(Solver& s, int d, V& reason);

  /* Generating reasons */
  Var eqi(Solver &s, int d) const;
  Var leqi(Solver &s, int d) const;

  // for when we know that d is in the original domain
  Var eqiUnsafe(Solver &s, int d) const;
  Var leqiUnsafe(Solver &s, int d) const;

  /* Generating reasons, convenience

     prefix r_ means 'reason'
     prefix e_ means 'effect'

     Note when constructing a reason, there can be at most one e_
     literal.
  */
  Lit r_geq(Solver &s, int d) const;
  Lit r_leq(Solver &s, int d) const;
  Lit r_neq(Solver &s, int d) const;
  Lit r_eq(Solver& s, int d) const;

  Lit r_min(Solver &s) const;
  Lit r_max(Solver &s) const;
  Lit r_eq(Solver &s) const;

  Lit e_geq(Solver &s, int d) const;
  Lit e_leq(Solver &s, int d) const;
  Lit e_neq(Solver &s, int d) const;
  Lit e_eq(Solver &s, int d) const;
};

bool operator==(cspvar x1, cspvar x2);
bool operator<(cspvar x1, cspvar x2);

// a stub for a constraint which goes in the priority queue.
struct consqueue
{
  int next, prev; // indices into consqs
  cons *c;
  int priority;

  consqueue() : next(-1), prev(-1), c(0L), priority(0) {}
  consqueue(cons *pc, int pp) : next(-1), prev(-1), c(pc), priority(pp) {}
};

/* A class for handling deferred explanations. Constraints subclass
   this and store any hints required to generate the actual
   explanation, or derive from both cons and explainer and use no
   other hints. In the former case, the explainer derived class will
   need to have a way to access the constraint that caused the
   pruning, for example by storing a pointer to it.

   Note that the constraint owns explainers. However, the solver will
   only use an explainer once for each literal (by calling explain())
   and the constraint is allowed to use this information to get rid of
   the explainer or reuse it. If an explainer's explain() is not
   called, release() will be called. At the time that the explainer is
   used as the reason, use() is called. In other words, #calls to
   use() = #calls to explain + #calls to release() and when that
   equilibrium is reached, the solver contains no more references to
   the explainer.
*/
class explainer
{
public:
    virtual ~explainer() {}

    /* The function actually called by the solver to generate the
       explanation. It makes c contain a clause that could be used as
       the reason for forcing p (this means that c includes p) */
    virtual void explain(Solver& s, Lit p, vec<Lit>& c) = 0;

    /* Called by the solver when it uses this explainer as the reason
       for pruning during propagation */
    virtual void use() = 0;
    /* Called by the solver when this explainer is no longer useful for
       a literal (meaning release may be called as many times as this
       object was used as the reason for a literal) */
    virtual void release() = 0;

    /* printer */
    virtual std::ostream& print(std::ostream& os)
    {
        return (os << "explainer at @" << this);
    }
};

inline std::ostream& operator<<(std::ostream& os, explainer& e)
{
    return e.print(os);
}

using explanation_ptr = pointer_variant<Clause, explainer>;

static Clause* const INVALID_CLAUSE{(Clause*)(0x1 << explanation_ptr::freebits)};
static Clause* const NO_REASON{nullptr};


/* The class from which all propagators derive */
class cons
{
  friend class Solver;
  int cqidx;
  int priority;
 public:
  cons() : cqidx(-1), priority(0) {}
  virtual ~cons() {}

  /* propagate, setting a clause or an explainer as reason for each
   * pruning. In case of failure, return a clause that describes the
   * reason.
   *
   * wake() is executed during propagation of a single variable, while
   * propagate() is a delayed call. wake_advised() is called when the
   * constraint gives a non-NULL piece of advice when it calls
   * Solver::wake_on_*
   */
  virtual Clause *wake_advised(Solver&, Lit, void *);
  virtual Clause *wake(Solver&, Lit);
  virtual Clause *propagate(Solver&);

  /* Clone this constraint and post it to othersolver. Useful for
   * debugging clauses */
  virtual void clone(Solver &othersolver);

  /* human readable representation of the constraint for debugging */
  virtual std::ostream& print(Solver& s, std::ostream& os) const;
  /* as above, but also any state (var domains, etc) */
  virtual std::ostream& printstate(Solver& s, std::ostream& os) const;

  /* Set scheduling priority. This has no effect if the propagator has
     been scheduled already, so it must be called in the constructor
     of a constraint, before it calls schedule_on_*
  */
  void set_priority(int p) {
    assert(0 <= p && p <= MAX_PRIORITY);
    priority = p;
  }

  /* Allows allocation optimizations similar to the one for Clause,
     which require that we call free(), not delete. */
  virtual void dispose() { delete this; }
};

typedef std::pair<cons*, void*> wake_stub;

// all the fixed or non-backtracked data of a csp variable
class cspvar_fixed
{
  friend class Solver;

  // dynamic data
  int min;
  int max;
  int dsize;

  // static data
  int omin; // min and max in the *original* domain
  int omax;
  Var firstbool;


  /* a cons may either wake immediately (like an ilog demon or a
     gecode advisor) when we process a literal, or it may be scheduled
     for execution later (although unlike gecode advisors, it may do
     propagation and cause failure in its wake()). When waking, it can
     also associate with the event a piece of advice to itself (e.g.,
     in n-ary propagators, the index of the variable that changed).
  */
  vec< wake_stub >
    wake_on_dom,
    wake_on_lb,
    wake_on_ub,
    wake_on_fix;
  vec<int>
    schedule_on_dom,
    schedule_on_lb,
    schedule_on_ub,
    schedule_on_fix;

  vec<Clause*> ps1, ps2, ps3, ps4;

  // accessing the propositional encoding
  bool ind(int i) const { return i >= omin && i <= omax; }
  Var eqi(int i) const { return ind(i)
      ? firstbool+2*(i-omin)
      : var_Undef; }
  Var leqi(int i) const { return ind(i)
      ? firstbool+2*(i-omin)+1
      : var_Undef; }

  Var eqiUnsafe(int i) const { return firstbool + 2*(i-omin); }
  Var leqiUnsafe(int i) const { return firstbool + 2*(i-omin)+1; }

  void copyclauses(vec<Clause*>& target, vec<Clause*> const& source)
  {
    target.growTo(source.size(), 0L);
    for(int i = 0; i != source.size(); ++i) {
      if( source[i] && source[i] != INVALID_CLAUSE )
        target[i] = Clause_new(*source[i]);
    }
  }
public:
  cspvar_fixed() {}
  cspvar_fixed(cspvar_fixed& f) :
    omin(f.omin), omax(f.omax), firstbool(f.firstbool)
  {
    copyclauses(ps1, f.ps1);
    copyclauses(ps2, f.ps2);
    copyclauses(ps3, f.ps3);
    copyclauses(ps4, f.ps4);
  }
};

struct domevent
{
  cspvar x;
  enum event_type { NEQ, EQ, GEQ, LEQ, NONE };
  event_type type;
  int d;

  domevent() :  x(-1), type(NONE), d(-1) {}
  domevent(cspvar px, event_type ptype, int pd) :
     x(px), type(ptype), d(pd) {}
};

inline bool noevent(domevent d) { return d.type == domevent::NONE; }
inline
const char *opstring(domevent::event_type t) {
  switch(t) {
  case domevent::EQ: return "==";
  case domevent::NEQ: return "!=";
  case domevent::LEQ: return "<=";
  case domevent::GEQ: return ">=";
  case domevent::NONE: return "#$%#^%";
  }
  return 0L; // gcc
}

//==================================================
// set vars

class setvar
{
  friend class Solver;
  int _id;
public:
  setvar() : _id(-1) {}
  explicit setvar(int id) : _id(id) {}

  int id() const { return _id; }
  bool valid() const { return _id >= 0; }

  // the min and max of its universe
  int umin(Solver &s) const;
  int umax(Solver &s) const;

  cspvar card(Solver &s) const;

  /* This is possibly misleading, but the names mean:
     includes d <=> d in lb
     excludes d <=> d notin ub
  */
  bool includes(Solver &s, int d) const;
  bool excludes(Solver &s, int d) const;

  /* Modify the domain */
  Clause *exclude(Solver &s, int d, Clause *c);
  Clause *exclude(Solver &s, int d, vec<Lit>& reason);
  Clause *exclude(Solver &s, int d, cons *c);
  Clause *excludef(Solver &s, int d, vec<Lit>& reason);

  Clause *include(Solver &s, int d, Clause *c);
  Clause *include(Solver &s, int d, vec<Lit>& reason);
  Clause *include(Solver &s, int d, cons *c);
  Clause *includef(Solver &s, int d, vec<Lit>& reason);

  /* Generate explanations */
  Var ini(Solver &s, int d) const;

  Lit r_ini(Solver &s, int d) const;
  Lit r_exi(Solver &s, int d) const;

  Lit e_ini(Solver &s, int d) const;
  Lit e_exi(Solver &s, int d) const;
};

bool operator==(setvar x1, setvar x2);

class setvar_data
{
  friend class Solver;

  int min;
  int max;
  Var firstbool;
  cspvar _card;

  vec< wake_stub > wake_on_in, wake_on_ex;
  vec<int> schedule_on_in, schedule_on_ex;

  bool inu(int i) const { return i >= min && i <= max; }
  Var ini(int i) const { return inu(i)
      ? firstbool + i - min
      : var_Undef; }
  cspvar card() const { return _card; }
public:
  setvar_data() {}
  setvar_data(setvar_data& s) :
    min(s.min), max(s.max), firstbool(s.firstbool),
    _card(s._card)
  {}
};

struct setevent
{
  setvar x;
  enum event_type { IN, EX, NONE };
  event_type type;
  int d;

  setevent() : x(0), type(NONE), d(0) {}
  setevent(setvar px, event_type ptype, int pd) :
    x(px), type(ptype), d(pd) {}
};

inline bool noevent(setevent d) { return d.type == setevent::NONE; }
inline
const char *opstring(setevent::event_type t) {
  switch(t) {
  case setevent::IN: return "in";
  case setevent::EX: return "notin";
  case setevent::NONE: return "%^%$#";
  }
  return 0L;
}

//==================================================
// output

struct cons_printer {
  Solver& _s;
  cons const & _c;
  cons_printer(Solver& s, cons const& c) : _s(s), _c(c) {}
};

inline
std::ostream& operator<<(std::ostream& os, cons_printer cp)
{
  return cp._c.print(cp._s, os);
}

struct cons_state_printer {
  Solver& _s;
  cons const & _c;
  cons_state_printer(Solver& s, cons const& c) : _s(s), _c(c) {}
};

inline
std::ostream& operator<<(std::ostream& os, cons_state_printer cp)
{
  return cp._c.printstate(cp._s, os);
}

struct var_printer {
  Solver& _s;
  Var _v;
  var_printer(Solver& s, Var v) : _s(s), _v(v) {}
};

std::ostream& operator<<(std::ostream& os, var_printer vp);

struct cspvar_printer {
  Solver& _s;
  cspvar _v;
  cspvar_printer(Solver& s, cspvar v) : _s(s), _v(v) {}
};

std::ostream& operator<<(std::ostream& os, cspvar_printer vp);

struct setvar_printer {
  Solver& _s;
  setvar _v;
  setvar_printer(Solver& s, setvar v) : _s(s), _v(v) {}
};

std::ostream& operator<<(std::ostream& os, setvar_printer vp);

struct lit_printer {
  Solver& _s;
  Lit _p;
  lit_printer(Solver& s, Lit p) : _s(s), _p(p) {}
};

std::ostream& operator<<(std::ostream& os, lit_printer lp);

struct domevent_printer {
  Solver& _s;
  domevent _p;
  domevent_printer(Solver& s, domevent p) : _s(s), _p(p) {}
};

std::ostream& operator<<(std::ostream& os, domevent_printer lp);

struct setevent_printer {
  Solver& _s;
  setevent _p;
  setevent_printer(Solver& s, setevent p) : _s(s), _p(p) {}
};

std::ostream& operator<<(std::ostream& os, setevent_printer lp);

inline bool is_invalid(Clause const* c) { return c == INVALID_CLAUSE; }
inline bool is_invalid(vec<Lit> const* v) { return false; }
inline bool is_invalid(std::vector<Lit> const* v) { return false; }

template<typename T>
struct clause_printer {
  Solver& _s;
  T const * _c;
  clause_printer(Solver& s, T const * c) : _s(s), _c(c) {}
};

template<typename T>
inline
std::ostream& operator<<(std::ostream& os, clause_printer<T> cp)
{
  if( !cp._c ) {
    os << "NO_REASON";
    return os;
  } else if (is_invalid(cp._c)) {
    os << "INVALID_CLAUSE";
    return os;
  }
  os << "( ";
  for(int i = 0; i != cp._c->size(); ++i)
    os << lit_printer(cp._s, (*cp._c)[i]) << " ";
  os << ")@" << cp._c;
  return os;
}

struct domain_as_range {
  Solver& _s;
  cspvar _x;
  domain_as_range(Solver &s, cspvar x) : _s(s), _x(x) {}
};

struct domain_as_set {
  Solver& _s;
  cspvar _x;
  domain_as_set(Solver &s, cspvar x) : _s(s), _x(x) {}
};

std::ostream& operator<<(std::ostream& os, domain_as_range x);
std::ostream& operator<<(std::ostream& os, domain_as_set x);

template<typename T>
inline
clause_printer<T> print(Solver& s, T const* c)
{
  return clause_printer<T>(s, c);
}

//==================================================
// exceptions

// thrown by the constructor of a constraint only
struct unsat : public std::exception
{
};

// thrown when we expect a var to be assigned
struct unassigned_var : public std::exception
{
  Solver const&_s;
  cspvar _x;

  unassigned_var(Solver const& s, cspvar x) : _s(s), _x(x) {}
  const char* what() const throw();
};

struct undefined_literal : public std::exception
{
  Solver const& _s;
  domevent _e;

  undefined_literal(Solver const& s, domevent e) :
    _s(s), _e(e) {}
  const char *what() const throw();
};

//==================================================
// convenience macros/inlines

#define DO_OR_THROW(action) do {                        \
    Clause *BOOST_PP_CAT(confl, __LINE__) = action;     \
    if( BOOST_PP_CAT(confl, __LINE__ ) ) throw unsat(); \
  } while(0)

#define DO_OR_RETURN(action) do {                       \
    Clause *BOOST_PP_CAT(confl, __LINE__ ) = action;    \
    if( BOOST_PP_CAT(confl, __LINE__ ) )                \
      return BOOST_PP_CAT(confl, __LINE__ );            \
  } while(0)

template<typename T>
inline
T* throw_if_null(T* t)
{
  if(!t) throw unsat();
  return t;
}

//======================================================================
// Provide a minimal vec<>-like interface to a (part of) a C
// array. The C array is meant to be on the stack, possibly allocated
// with alloca. This is a very dangerous class: the user must ensure
// that the array passed is large enough to accept a push(), if it is
// used.

template<typename T>
class stack_array
{
  size_t _n;
  T* _data;
public:
  stack_array() : _n(0), _data(0L) {}
  stack_array(size_t n, T* data) : _n(n), _data(data) {}

  T& operator[](size_t i) { return _data[i]; }
  T const& operator[](size_t i) const { return _data[i]; }
  // hmph, this should return size_t but vec<> returns int
  int size() const { return _n; }

  void push(T const& t) { _data[_n] = t; ++_n; }
  void pop() { --_n; }

  void shrink(int i) { _n -= i; }
};

} // namespace minicsp

#endif
