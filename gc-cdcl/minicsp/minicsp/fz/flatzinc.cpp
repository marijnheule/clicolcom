/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*****************************************************************************
  minicsp

Copyright 2010 George Katsirelos

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

Parts of this file were distributed under the following license as part
of Gecode.

  Main authors:
     Guido Tack <tack@gecode.org>

  Copyright:
     Guido Tack, 2007

  Last modified:
     $Date: 2010-05-11 12:33:38 +0200 (Tue, 11 May 2010) $ by $Author: tack $
     $Revision: 10940 $

  This file is part of Gecode, the generic constraint
  development environment:
     http://www.gecode.org

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*****************************************************************************/


#include "flatzinc.hpp"
#include "registry.hpp"
#include <iomanip>

#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;

namespace minicsp {

namespace FlatZinc {
  inline
  set<int> setrange(int min, int max)
  {
    set<int> rv;
    for(int i = min; i <= max; ++i)
      rv.insert(i);
    return rv;
  }

  set<int> vs2is(IntVarSpec* vs) {
    if (vs->assigned) {
      return setrange(vs->i,vs->i);
    }
    if (vs->domain()) {
      AST::SetLit* sl = vs->domain.some();
      if (sl->interval) {
        return setrange(sl->min, sl->max);
      } else {
        set<int> rv;
        for (int i=sl->s.size(); i--;)
          rv.insert(sl->s[i]);
        return rv;
      }
    }
    return setrange(-1000, 1000);
  }

  int vs2bsl(BoolVarSpec* bs) {
    if (bs->assigned) {
      return bs->i;
    }
    if (bs->domain()) {
      AST::SetLit* sl = bs->domain.some();
      assert(sl->interval);
      return std::min(1, std::max(0, sl->min));
    }
    return 0;
  }

  int vs2bsh(BoolVarSpec* bs) {
    if (bs->assigned) {
      return bs->i;
    }
    if (bs->domain()) {
      AST::SetLit* sl = bs->domain.some();
      assert(sl->interval);
      return std::max(0, std::min(1, sl->max));
    }
    return 1;
  }

  int ann2ivarsel(AST::Node* ann) {
    if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
      // if (s->id == "input_order")
      //   return TieBreakVarBranch<IntVarBranch>(INT_VAR_NONE);
    }
    std::cerr << "Warning, ignored search annotation: ";
    ann->print(std::cerr);
    std::cerr << std::endl;
    return 0;
  }

  int ann2ivalsel(AST::Node* ann) {
    if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
      // if (s->id == "indomain_min")
      //   return INT_VAL_MIN;
    }
    std::cerr << "Warning, ignored search annotation: ";
    ann->print(std::cerr);
    std::cerr << std::endl;
    return 0;
  }

  int ann2asnivalsel(AST::Node* ann) {
    if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
      // if (s->id == "indomain_min")
      //   return INT_ASSIGN_MIN;
    }
    std::cerr << "Warning, ignored search annotation: ";
    ann->print(std::cerr);
    std::cerr << std::endl;
    return 0;
  }


  FlatZincModel::FlatZincModel(Solver &s)
    : solver(s),
      intVarCount(-1), boolVarCount(-1), setVarCount(-1), _optVar(-1),
      _solveAnnotations(NULL),
      findall(false)
  {}

  void
  FlatZincModel::init(int intVars, int boolVars, int setVars) {
    intVarCount = 0;
    iv = IntVarArray(intVars);
    iv_introduced = std::vector<bool>(intVars);
    iv_boolalias = std::vector<int>(intVars);
    boolVarCount = 0;
    bv = BoolVarArray(boolVars);
    bv_introduced = std::vector<bool>(boolVars);
    setVarCount = 0;
    sv = SetVarArray(setVars);
    sv_introduced = std::vector<bool>(setVars);
  }

  void
  FlatZincModel::newIntVar(IntVarSpec* vs) {
    if (vs->alias) {
      iv[intVarCount++] = iv[vs->i];
    } else {
      set<int> domain = vs2is(vs);
      cspvar x = solver.newCSPVar(*domain.begin(), *domain.rbegin());
      iv[intVarCount++] = x;
      int prev = *domain.begin();
      for(set<int>::const_iterator i = domain.begin(), end = domain.end();
          i != end; ++i) {
        if( *i > prev+1 ) {
          for(int q = prev+1; q != *i; ++q)
            x.remove(solver, q, NO_REASON);
        }
        prev = *i;
      }
    }
    iv_introduced[intVarCount-1] = vs->introduced;
    iv_boolalias[intVarCount-1] = -1;
  }

  void
  FlatZincModel::newSetVar(SetVarSpec* vs) {
    if (vs->alias) {
      sv[intVarCount++] = sv[vs->i];
    } else if( vs->assigned) {
      assert(vs->upperBound());
      AST::SetLit* vsv = vs->upperBound.some();
      if (vsv->interval) {
        setvar x = solver.newSetVar(vsv->min, vsv->max);
        sv[setVarCount++] = x;
        for(int i = vsv->min; i <= vsv->max; ++i)
          x.include(solver, i, NO_REASON);
      } else {
        if( vsv->s.empty() ) {
          setvar x = solver.newSetVar( 0, 0 );
          x.exclude(solver, 0, NO_REASON);
          sv[setVarCount++] = x;
        } else {
          int umin = vsv->s[0], umax = vsv->s[0];
          for(size_t i = 1; i != vsv->s.size(); ++i) {
            umin = std::min(umin, vsv->s[i]);
            umax = std::max(umax, vsv->s[i]);
          }
          setvar x = solver.newSetVar(umin, umax);
          sv[setVarCount++] = x;
          for(size_t i = 0; i != vsv->s.size(); ++i)
            x.include(solver, vsv->s[i], NO_REASON);
          for(int i = x.umin(solver), iend = x.umax(solver); i <= iend; ++i)
            if( !x.includes(solver, i) )
              x.exclude(solver, i, NO_REASON);
        }
      }
    } else if( vs->upperBound() ) {
      AST::SetLit* vsv = vs->upperBound.some();
      setvar x = solver.newSetVar(vsv->min, vsv->max);
      sv[setVarCount++] = x;
      if( !vsv->interval ) {
        int prev = vsv->min;
        for(size_t i = 0; i != vsv->s.size(); ++i) {
          if( vsv->s[i] > prev+1 ) {
            for(int q = prev+1; q != vsv->s[i]; ++q)
              x.exclude(solver, q, NO_REASON);
          }
          prev = vsv->s[i];
        }
      } // otherwise everything is unset and we are done here
    } else {
      // completely free
      setvar x = solver.newSetVar(-1000, 1000);
      sv[setVarCount++] = x;
    }
    sv_introduced[setVarCount-1] = vs->introduced;
  }

  void
  FlatZincModel::aliasBool2Int(int iv, int bv) {
    iv_boolalias[iv] = bv;
  }
  int
  FlatZincModel::aliasBool2Int(int iv) {
    return iv_boolalias[iv];
  }

  void
  FlatZincModel::newBoolVar(BoolVarSpec* vs) {
    if (vs->alias) {
      bv[boolVarCount++] = bv[vs->i];
    } else {
      bv[boolVarCount++] = solver.newCSPVar(vs2bsl(vs), vs2bsh(vs));
    }
    bv_introduced[boolVarCount-1] = vs->introduced;
  }

  void
  FlatZincModel::postConstraint(const ConExpr& ce, AST::Node* ann) {
    try {
      registry().post(solver, *this, ce, ann);
    } catch (AST::TypeError& e) {
      throw FlatZinc::Error("Type error", e.what());
    }
  }

  void flattenAnnotations(AST::Array* ann, std::vector<AST::Node*>& out) {
      for (unsigned int i=0; i<ann->a.size(); i++) {
        if (ann->a[i]->isCall("seq_search")) {
          AST::Call* c = ann->a[i]->getCall();
          if (c->args->isArray())
            flattenAnnotations(c->args->getArray(), out);
          else
            out.push_back(c->args);
        } else {
          out.push_back(ann->a[i]);
        }
      }
  }

  void
  FlatZincModel::createBranchers(AST::Node* ann, bool ignoreUnknown,
                                 std::ostream& err) {
    if (ann) {
      err << "Warning, ignored search annotation: ";
      ann->print(err);
      err << std::endl;
    }
  }

  AST::Array*
  FlatZincModel::solveAnnotations(void) const {
    return _solveAnnotations;
  }

  void
  FlatZincModel::solve(AST::Array* ann) {
    _method = SAT;
    _solveAnnotations = ann;
  }

  void
  FlatZincModel::minimize(int var, AST::Array* ann) {
    _method = MIN;
    _optVar = var;
    _solveAnnotations = ann;
    // Branch on optimization variable to ensure that it is given a value.
    AST::Array* args = new AST::Array(4);
    args->a[0] = new AST::Array(new AST::IntVar(_optVar));
    args->a[1] = new AST::Atom("input_order");
    args->a[2] = new AST::Atom("indomain_min");
    args->a[3] = new AST::Atom("complete");
    AST::Call* c = new AST::Call("int_search", args);
    if (!ann)
      ann = new AST::Array(c);
    else
      ann->a.push_back(c);
  }

  void
  FlatZincModel::maximize(int var, AST::Array* ann) {
    _method = MAX;
    _optVar = var;
    _solveAnnotations = ann;
    // Branch on optimization variable to ensure that it is given a value.
    AST::Array* args = new AST::Array(4);
    args->a[0] = new AST::Array(new AST::IntVar(_optVar));
    args->a[1] = new AST::Atom("input_order");
    args->a[2] = new AST::Atom("indomain_min");
    args->a[3] = new AST::Atom("complete");
    AST::Call* c = new AST::Call("int_search", args);
    if (!ann)
      ann = new AST::Array(c);
    else
      ann->a.push_back(c);
  }

  FlatZincModel::~FlatZincModel(void) {
    delete _solveAnnotations;
  }

  void
  FlatZincModel::run(std::ostream& out, const Printer& p) {
    using std::setw;
    using std::setfill;
    // solve it
    bool sat = false, next;
    do {
      next = solver.solve();
      sat = sat || next;
      if( next ) {
        print(out, p);
        out << setw(10) << setfill('-') << "-" << "\n";
        switch (_method) {
        case MIN:
        case MAX:
          try {
            constrain();
          } catch( unsat ) {
            // we have already derived a lower bound equal to this
            // solution
            next = false;
          }
          break;
        case SAT:
          next = findall;
          if( findall ) {
            try {
              solver.excludeLast();
            } catch( unsat ) {
              // no more solutions :( Poor us
              next = false;
            }
          }
          break;
        }
      }
    } while(next);
    if( sat )
      out << setw(10) << setfill('=') << '=' << "\n";
    else
      out << setw(5) << setfill('=') << '='
          << "UNSATISFIABLE" << setw(5) << '=' << "\n";
  }

  void
  FlatZincModel::constrain() {
    int opt = solver.cspModelValue(iv[_optVar]);
    if (_method == MIN) {
      iv[_optVar].setmax( solver, opt-1, NO_REASON );
    } else if (_method == MAX) {
      iv[_optVar].setmin( solver, opt+1, NO_REASON );
    }
  }

  FlatZincModel::Meth
  FlatZincModel::method(void) const {
    return _method;
  }

  int
  FlatZincModel::optVar(void) const {
    return _optVar;
  }

  void
  FlatZincModel::print(std::ostream& out, const Printer& p) const {
    p.print(out, solver, iv, bv, sv);
  }

  void
  Printer::init(AST::Array* output) {
    _output = output;
  }

  void
  Printer::printElem(std::ostream& out,
                     Solver& solver,
                     AST::Node* ai,
                     const IntVarArray& iv,
                     const BoolVarArray& bv,
                     const SetVarArray& sv
                       ) const {
    int k;
    if (ai->isInt(k)) {
      out << k;
    } else if (ai->isIntVar()) {
      pair<int, int> v = solver.cspModelRange(iv[ai->getIntVar()]);
      if( v.first == v.second )
        out << v.first;
      else
        out << v.first << ".." << v.second;
    } else if (ai->isBoolVar()) {
      pair<int, int> v = solver.cspModelRange(bv[ai->getBoolVar()]);
      if (v.first == 1) {
        out << "true";
      } else if (v.second == 0) {
        out << "false";
      } else {
        out << "false..true";
      }
    } else if( ai->isSetVar()) {
      setvar x = sv[ai->getSetVar()];
      pair< set<int>, set<int> > const& v = solver.cspSetModel(x);
      set<int> const& lb = v.first;
      set<int> const& ub = v.second;
      assert( lb == ub );
      out << "{";
      for( set<int>::const_iterator i = ub.begin(); i != ub.end(); ++i) {
        if( i != ub.begin() ) out << ", ";
        out << *i;
      }
      out << "}";
    } else if (ai->isBool()) {
      out << (ai->getBool() ? "true" : "false");
    } else if (ai->isSet()) {
      AST::SetLit* s = ai->getSet();
      if (s->interval) {
        out << s->min << ".." << s->max;
      } else {
        out << "{";
        for (unsigned int i=0; i<s->s.size(); i++) {
          out << s->s[i] << (i < s->s.size()-1 ? ", " : "}");
        }
      }
    } else if (ai->isString()) {
      std::string s = ai->getString();
      for (unsigned int i=0; i<s.size(); i++) {
        if (s[i] == '\\' && i<s.size()-1) {
          switch (s[i+1]) {
          case 'n': out << "\n"; break;
          case '\\': out << "\\"; break;
          case 't': out << "\t"; break;
          default: out << "\\" << s[i+1];
          }
          i++;
        } else {
          out << s[i];
        }
      }
    }
  }

  void
  Printer::print(std::ostream& out,
                 Solver& solver,
                 const IntVarArray& iv,
                 const BoolVarArray& bv,
                 const SetVarArray& sv) const {
    if (_output == NULL)
      return;
    for (unsigned int i=0; i< _output->a.size(); i++) {
      AST::Node* ai = _output->a[i];
      if (ai->isArray()) {
        AST::Array* aia = ai->getArray();
        int size = aia->a.size();
        out << "[";
        for (int j=0; j<size; j++) {
          printElem(out,solver, aia->a[j],iv,bv,sv);
          if (j<size-1)
            out << ", ";
        }
        out << "]";
      } else {
        printElem(out,solver,ai,iv,bv,sv);
      }
    }
  }

  Printer::~Printer(void) {
    delete _output;
  }

}

} //namespace minicsp

// STATISTICS: flatzinc-any
