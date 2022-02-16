/******************************************************************************
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

Parts of this file were distributed under the following license as
part of the xcsp distribution

 Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
******************************************************************************/

#ifndef _MINICSP_CALLBACK_H_
#define _MINICSP_CALLBACK_H_

#include <vector>
#include <map>
#include "CSPParserCallback.hh"
#include "minicsp/core/cons.hpp"

using namespace minicsp;

struct unsupported {};

struct minicsp_rel {
  RelType type;
  std::vector< std::vector<int> > tuples;
};

struct minicsp_pred {
  std::map<std::string, int> par;
  CSPXMLParser::AST *tree;
};

namespace CSPXMLParser
{
  using namespace std;

  class FixedValuation : public VariableValuation {
    vector<int> _vals;
  public:
    FixedValuation(vector<int> const& vals) : _vals(vals) {}
    FunctionValue getVarValue(unsigned id) const { return _vals[id]; }
  };

// definition of the functions which are called by the parser when it
// reads new information in the file. These functions are in charge of
// creating the data structure that the solver needs to do his job
class minicsp_callback : public CSPParserCallback
{
  Solver &_solver;

  typedef vector< pair<int, int> > domain_t;
  map<string, domain_t> domains;
  string current_domain_name;
  domain_t current_domain;

  map<string, cspvar> vars;
  map<int, cspvar> constants;

  map<string, minicsp_rel> relations;
  string current_relation_name;
  minicsp_rel current_relation;

  map<string, minicsp_pred> preds;
  string current_pred_name;
  minicsp_pred current_pred;

  string constraintReference;
  vector<cspvar> current_scope;


  vector<cspvar> ast2vararray(AST const& ast)
  {
    return ast2vararray( dynamic_cast<ASTList const&>(*this) );
  }

  cspvar ast2var(AST const& ast)
  {
    if( ast.isVar() )
      return vars[ ast.getVarName() ];
    else {
      int c = ast.getInteger();
      if( constants.find(c) == constants.end() ) {
        constants[c] = _solver.newCSPVar(c, c);
      }
      return constants[c];
    }
  }

  vector<cspvar> ast2vararray(ASTList const& ast)
  {
    vector<cspvar> v(ast.size());
    for(int i = 0; i != ast.size(); ++i) {
      v[i] = ast2var(ast[i]);
    }
    return v;
  }

  void post_table(string const& reference,
                  ASTList const& args)
  {
    minicsp_rel const& rel = relations[reference];
    if( rel.type == REL_SUPPORT )
      post_positive_table(_solver, current_scope, rel.tuples);
    else
      post_negative_table(_solver, current_scope, rel.tuples);
  }

  void post_expression_table(AST *tree, vector<cspvar> const& vars,
                             size_t idx,
                             vector<int> &curtuple,
                             vector< vector<int> > & curtable)
  {
    if( idx == vars.size() ) {
      if(tree->value(FixedValuation(curtuple)))
        curtable.push_back(curtuple);
      return;
    }
    for(int q = vars[idx].min(_solver); q <= vars[idx].max(_solver); ++q) {
      if( vars[idx].indomain(_solver, q) ) {
        curtuple[idx] = q;
        post_expression_table(tree, vars, idx+1, curtuple, curtable);
      }
    }
    if( idx == 0 )
      post_positive_table(_solver, vars, curtable);
  }

  cspvar post_expression(C_AST *ctree, vector<cspvar> const& vars, bool root)
  {
    cspvar rv;
    switch( ctree->type ) {
    case VAR: {
      assert(!root);
      C_AST_VarNode *vn = (C_AST_VarNode*)ctree;
      rv = vars[vn->idVar];
      break;
    }
    case F_EQ: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_eq(_solver, x1, x2, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_eq_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_NE: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_neq(_solver, x1, x2, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_neq_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_GE: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_leq(_solver, x2, x1, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_geq_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_GT: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_less(_solver, x2, x1, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_gt_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_LE: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_leq(_solver, x1, x2, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_leq_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_LT: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar x1 = post_expression(fn->args[0], vars, false);
      cspvar x2 = post_expression(fn->args[1], vars, false);
      if( root )
        post_less(_solver, x1, x2, 0);
      else {
        rv = _solver.newCSPVar(0,1);
        post_less_re(_solver, x1, x2, 0, rv);
      }
      break;
    }
    case F_OR: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( root ) {
        vec<Lit> ps;
        for(int i = 0; i != fn->nbarg; ++i ) {
          cspvar arg = post_expression(fn->args[i], vars, false);
          ps.push( arg.r_eq(_solver, 0) );
        }
        _solver.addClause(ps);
      } else {
        vec<Lit> ps;
        rv = _solver.newCSPVar(0,1);
        ps.push( rv.e_eq(_solver, 0) );
        for(int i = 0; i != fn->nbarg; ++i) {
          cspvar arg = post_expression(fn->args[i], vars, false);

          vec<Lit> ps1;
          ps1.push( rv.r_eq(_solver, 0) );
          ps1.push( arg.e_eq(_solver, 0) );
          _solver.addClause(ps1);

          ps.push( arg.r_eq(_solver, 0) );
        }
        _solver.addClause(ps);
      }
      break;
    }
    case F_AND: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( root ) {
        for(int i = 0; i != fn->nbarg; ++i )
          post_expression(fn->args[i], vars, true);
      } else {
        vec<Lit> ps;
        rv = _solver.newCSPVar(0,1);
        ps.push( rv.r_eq(_solver, 0) );
        for(int i = 0; i != fn->nbarg; ++i) {
          cspvar arg = post_expression(fn->args[i], vars, false);

          vec<Lit> ps1;
          ps1.push( rv.e_eq(_solver, 0) );
          ps1.push( arg.r_eq(_solver, 0) );
          _solver.addClause(ps1);

          ps.push( arg.e_eq(_solver, 0) );
        }
        _solver.addClause(ps);
      }
      break;
    }
    case F_NOT: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      assert(fn->nbarg == 1);
      if( root ) {
        cspvar x = post_expression(fn->args[0], vars, true);
        DO_OR_THROW(x.assign(_solver, 0, NO_REASON));
      } else {
        vec<Lit> ps;
        rv = _solver.newCSPVar(0,1);
        cspvar arg = post_expression(fn->args[0], vars, false);
        vec<Lit> ps1, ps2;
        ps1.push( rv.r_eq(_solver, 0) );
        ps1.push( arg.e_neq(_solver, 0) );
        ps2.push( rv.r_neq(_solver, 0) );
        ps2.push( arg.e_eq(_solver, 0) );
        _solver.addClause(ps1);
        _solver.addClause(ps2);
      }
      break;
    }
    case F_IFF: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      assert(fn->nbarg == 2);
      if( root ) {
        cspvar arg1 = post_expression(fn->args[0], vars, false),
          arg2 = post_expression(fn->args[1], vars, false);
        vec<Lit> ps1, ps2;
        ps1.push( arg1.r_eq(_solver, 0) );
        ps1.push( arg2.e_eq(_solver, 0) );
        ps2.push( arg1.r_neq(_solver, 0) );
        ps2.push( arg2.e_neq(_solver, 0) );
        _solver.addClause(ps1);
        _solver.addClause(ps2);
      } else {
        vec<Lit> ps;
        rv = _solver.newCSPVar(0,1);
        cspvar arg1 = post_expression(fn->args[0], vars, false),
          arg2 = post_expression(fn->args[1], vars, false);
        vec<Lit> ps1, ps2, ps3, ps4;
        ps1.push( rv.r_neq(_solver, 0) );
        ps1.push( arg1.r_eq(_solver, 0) );
        ps1.push( arg2.e_eq(_solver, 0) );
        ps2.push( rv.r_neq(_solver, 0) );
        ps2.push( arg1.r_neq(_solver, 0) );
        ps2.push( arg2.e_neq(_solver, 0) );

        ps3.push( rv.r_eq(_solver, 0) );
        ps3.push( arg1.r_eq(_solver, 0) );
        ps3.push( arg2.e_neq(_solver, 0) );
        ps4.push( rv.r_eq(_solver, 0) );
        ps4.push( arg1.r_neq(_solver, 0) );
        ps4.push( arg2.e_eq(_solver, 0) );
        _solver.addClause(ps1);
        _solver.addClause(ps2);
        _solver.addClause(ps3);
        _solver.addClause(ps4);
      }
      break;
    }

    case F_XOR: {
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      assert(fn->nbarg == 2);
      if( root ) {
        cspvar arg1 = post_expression(fn->args[0], vars, false),
          arg2 = post_expression(fn->args[1], vars, false);
        vec<Lit> ps1, ps2;
        ps1.push( arg1.r_eq(_solver, 0) );
        ps1.push( arg2.e_neq(_solver, 0) );
        ps2.push( arg1.r_neq(_solver, 0) );
        ps2.push( arg2.e_eq(_solver, 0) );
        _solver.addClause(ps1);
        _solver.addClause(ps2);
      } else {
        vec<Lit> ps;
        rv = _solver.newCSPVar(0,1);
        cspvar arg1 = post_expression(fn->args[0], vars, false),
          arg2 = post_expression(fn->args[1], vars, false);
        vec<Lit> ps1, ps2, ps3, ps4;
        ps1.push( rv.r_neq(_solver, 0) );
        ps1.push( arg1.r_eq(_solver, 0) );
        ps1.push( arg2.e_neq(_solver, 0) );
        ps2.push( rv.r_neq(_solver, 0) );
        ps2.push( arg1.r_neq(_solver, 0) );
        ps2.push( arg2.e_eq(_solver, 0) );

        ps3.push( rv.r_eq(_solver, 0) );
        ps3.push( arg1.r_eq(_solver, 0) );
        ps3.push( arg2.e_eq(_solver, 0) );
        ps4.push( rv.r_eq(_solver, 0) );
        ps4.push( arg1.r_neq(_solver, 0) );
        ps4.push( arg2.e_neq(_solver, 0) );
        _solver.addClause(ps1);
        _solver.addClause(ps2);
        _solver.addClause(ps3);
        _solver.addClause(ps4);
      }
      break;
    }

    // function stuff
    case F_NEG: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 1 )
        throw unsupported();
      cspvar arg = post_expression(fn->args[0], vars, false);
      rv = _solver.newCSPVar(-arg.max(_solver), -arg.min(_solver));
      post_neg(_solver, arg, rv, 0);
      break;
    }
    case F_ABS: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 1 )
        throw unsupported();
      cspvar arg = post_expression(fn->args[0], vars, false);
      rv = _solver.newCSPVar(0, max(abs(arg.min(_solver)), abs(arg.max(_solver))));
      post_abs(_solver, arg, rv, 0);
      break;
    }
    case F_SUB: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      if( fn->nbarg != 2 )
        throw unsupported();
      cspvar arg1 = post_expression(fn->args[0], vars, false);
      cspvar arg2 = post_expression(fn->args[1], vars, false);
      int min = arg1.min(_solver) - arg2.max(_solver);
      int max = arg1.max(_solver) - arg2.min(_solver);
      rv = _solver.newCSPVar(min, max);
      vector<int> w(3);
      vector<cspvar> v(3);
      w[0] = 1;   v[0] = rv;
      w[1] = -1;  v[1] = arg1;
      w[2] = 1;   v[2] = arg2;
      post_lin_eq(_solver, v, w, 0);
      break;
    }
    case F_ADD: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      int min = 0, max = 0;
      vector<int> w;
      vector<cspvar> v;
      for(int q = 0; q != fn->nbarg; ++q) {
        cspvar arg = post_expression(fn->args[q], vars, false);
        w.push_back(-1);
        v.push_back( arg );
        min += arg.min(_solver);
        max += arg.max(_solver);
      }
      rv = _solver.newCSPVar(min, max);
      v.push_back(rv);
      w.push_back(1);
      post_lin_eq(_solver, v, w, 0);
      break;
    }
    case F_MUL: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      cspvar arg0 = post_expression(fn->args[0], vars, false);
      for(int q = 1; q != fn->nbarg; ++q) {
        cspvar arg1 = post_expression(fn->args[q], vars, false);
        int minv = min( min( arg0.min(_solver)*arg1.min(_solver),
                             arg0.min(_solver)*arg1.max(_solver)),
                        min( arg0.max(_solver)*arg1.min(_solver),
                             arg0.max(_solver)*arg1.max(_solver))),
          maxv = max( max( arg0.min(_solver)*arg1.min(_solver),
                           arg0.min(_solver)*arg1.max(_solver)),
                      max( arg0.max(_solver)*arg1.min(_solver),
                           arg0.max(_solver)*arg1.max(_solver)));
        cspvar res = _solver.newCSPVar(minv, maxv);
        post_mult(_solver, res, arg0, arg1);
        arg0 = res;
      }
      rv = arg0;
      break;
    }

    case F_MIN: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      cspvar arg0 = post_expression(fn->args[0], vars, false);
      for(int q = 1; q != fn->nbarg; ++q) {
        cspvar arg1 = post_expression(fn->args[q], vars, false);
        int minv = min( arg0.min(_solver), arg1.min(_solver) ),
          maxv = min(arg0.max(_solver), arg1.max(_solver));
        cspvar res = _solver.newCSPVar(minv, maxv);
        post_min(_solver, res, arg0, arg1);
        arg0 = res;
      }
      rv = arg0;
      break;
    }

    case F_MAX: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      cspvar arg0 = post_expression(fn->args[0], vars, false);
      for(int q = 1; q != fn->nbarg; ++q) {
        cspvar arg1 = post_expression(fn->args[q], vars, false);
        int minv = max( arg0.min(_solver), arg1.min(_solver) ),
          maxv = max(arg0.max(_solver), arg1.max(_solver));
        cspvar res = _solver.newCSPVar(minv, maxv);
        post_max(_solver, res, arg0, arg1);
        arg0 = res;
      }
      rv = arg0;
      break;
    }

    case F_IF: {
      assert(!root);
      C_AST_FxNode *fn = (C_AST_FxNode*)(ctree);
      assert( fn->nbarg == 3 );
      cspvar argif = post_expression(fn->args[0], vars, false);
      cspvar arg1 = post_expression(fn->args[1], vars, false);
      cspvar arg2 = post_expression(fn->args[2], vars, false);
      int minv = min(arg1.min(_solver), arg2.min(_solver)),
        maxv = max(arg1.max(_solver), arg2.max(_solver));
      rv = _solver.newCSPVar(minv, maxv);
      post_eq_re(_solver, rv, arg1, 0, argif.e_eq(_solver, 1));
      post_eq_re(_solver, rv, arg2, 0, argif.e_neq(_solver, 1));
      break;
    }

      // all the unimplemented stuff goes here
    case F_POW:
    case F_MOD:
    case F_DIV:
    case CST_BOOL:
    case CST_INT:
    case LIST:
    case DICT:
    default:
      throw unsupported();
    }
    return rv;
  }

  void post_expression(string const& reference,
                       ASTList const& args)
  {
    bool table_expr = false;
    if( current_scope.size() <= 2 )
      table_expr = true;
    if( !table_expr ) {
      C_AST *ctree = preds[reference].tree->makeCTree();
      post_expression(ctree, ast2vararray(args), true );
    } else {
      vector<cspvar> X = ast2vararray(args);
      vector<int> curtuple(X.size());
      vector< vector<int> > table;
      post_expression_table(preds[reference].tree, X,
                            0, curtuple, table);
    }
  }

  void post_alldiff(string const& reference,
                    ASTList const& args)
  {
    ::post_alldiff(_solver, current_scope);
  }

  void post_wsum(string const& reference,
                 ASTList const& args)
  {
    const AST &sum=args[0];
    const ASTSymb &op=dynamic_cast<ASTSymb const&>(args[1]);
    const AST &rhs=args[2];

    vector<cspvar> const& vars = current_scope;

    vector<int> w;
    vector<cspvar> v;

    for(int i=0;i<sum.size();++i) {
      w.push_back( sum[i]["coef"].getInteger() );
      v.push_back( vars[ sum[i]["var"].getVarId() ] );
    }

    int c = 0;

    if( rhs.isVar() ) {
      w.push_back(-1);
      v.push_back( vars[ rhs.getVarId() ] );
    } else {
      c = -rhs.getInteger();
    }

    switch(op.getType()) {
    case SYMB_EQ:
      post_lin_eq(_solver, v, w, c);
      break;
    case SYMB_NE:
      post_lin_neq(_solver, v, w, c);
      break;
    case SYMB_GE:
      for(int i = 0; i != sum.size(); ++i) w[i] = -w[i];
      post_lin_leq(_solver, v, w, c);
      break;
    case SYMB_GT:
      for(int i = 0; i != sum.size(); ++i) w[i] = -w[i];
      post_lin_less(_solver, v, w, c);
      break;
    case SYMB_LE:
      post_lin_leq(_solver, v, w, c);
      break;
    case SYMB_LT:
      post_lin_less(_solver, v, w, c);
      break;
    default:
      throw unsupported();
    }
  }

  void post_cumulative(string const& reference,
                       ASTList const& args)
  {
    const AST &tasks=args[0];
    const AST &limit=args[1];

    cout << "cumulative global constraint:\n  "
         << tasks.size() << " tasks are defined\n";

    for(int i=0;i<tasks.size();++i)
      {
        const AST &desc=tasks[i];
        cout << "    task " << i << "= ( ";
        if (desc.hasKey("origin"))
          {
            cout << "origin=";
            if (desc["origin"].isVar())
              cout << desc["origin"].getVarName();
            else
              if (desc["origin"].isInteger())
                cout << desc["origin"].getInteger();
              else
                cout << "?";

            cout << ' ';
          }

        if (desc.hasKey("duration"))
          {
            cout << "duration=";
            if (desc["duration"].isVar())
              cout << desc["duration"].getVarName();
            else
              if (desc["duration"].isInteger())
                cout << desc["duration"].getInteger();
              else
                cout << "?";

            cout << ' ';
          }

        if (desc.hasKey("end"))
          {
            cout << "end=";
            if (desc["end"].isVar())
              cout << desc["end"].getVarName();
            else
              if (desc["end"].isInteger())
                cout << desc["end"].getInteger();
              else
                cout << "?";

            cout << ' ';
          }

        if (desc.hasKey("height"))
          {
            cout << "height=";
            if (desc["height"].isVar())
              cout << desc["height"].getVarName();
            else
              if (desc["height"].isInteger())
                cout << desc["height"].getInteger();
              else
                cout << "?";

            cout << ' ';
          }

        cout << ")\n";
      }

    cout << "  limit=";
    if (limit.isVar())
      cout << limit.getVarName();
    else
      if (limit.isInteger())
        cout << limit.getInteger();
      else
        cout << "?";

    cout << endl;
  }

  void post_element(string const& reference,
                    ASTList const& args)
  {
    const AST &indexast=args[0];
    const AST &arrayast=args[1];
    const AST &valueast=args[2];

    cspvar index = ast2var(indexast);
    cspvar result = ast2var(valueast);
    vector<cspvar> array = ast2vararray(arrayast);
    // FIXME: find out if it is really 1-based indexing, as in the
    // catalog
    ::post_element(_solver, result, index, array, 1);
  }
public:
  minicsp_callback(Solver& s) : _solver(s) {
    posters["global:alldifferent"] = &minicsp_callback::post_alldiff;
    // to be uncommented when there is an actual implementation there
    posters["global:element"] = &minicsp_callback::post_element;
    posters["global:weightedsum"] = &minicsp_callback::post_wsum;
    //posters["global:cumulative"] = &minicsp_callback::post_cumulative;
  }

  virtual void beginDomain(const string & name, int idDomain, int nbValue)
  {
    current_domain_name = name;
  }

  void addDomainValue(int v)
  {
    current_domain.push_back( make_pair(v, v) );
  }

  virtual void addDomainValue(int first,int last)
  {
    current_domain.push_back( make_pair(first, last) );
  }

  virtual void endDomain()
  {
    domains[current_domain_name] = current_domain;
    current_domain.clear();
  }

  virtual void addVariable(const string & name, int idVar,
                           const string & domain, int idDomain)
  {
    domain_t const& d = domains[domain];
    cspvar x = _solver.newCSPVar(d.front().first, d.back().second);
    for(size_t i = 1; i != d.size(); ++i) {
      for(int j = d[i-1].second+1; j < d[i].first; ++j)
        x.remove(_solver, j, NO_REASON);
    }
    _solver.setCSPVarName(x, name);
    vars[name] = x;
  }

  virtual void beginRelation(const string & name, int idRel,
                             int arity, int nbTuples, RelType relType)
  {
    current_relation_name = name;
    current_relation.tuples.clear();
    current_relation.type = relType;

    if (relType != REL_SUPPORT && relType != REL_CONFLICT )
      throw unsupported();
  }

  virtual void addRelationTuple(int arity, int tuple[])
  {
    vector<int> t(tuple, tuple+arity);
    current_relation.tuples.push_back(t);
  }

  virtual void endRelation()
  {
    relations[current_relation_name] = current_relation;
    posters[current_relation_name] = &minicsp_callback::post_table;
  }

  virtual void beginPredicate(const string & name, int idPred)
  {
    current_pred_name = name;
    current_pred.par.clear();
  }

  virtual void addFormalParameter(int pos, const string & name,
                                  const string & type)
  {
    if( type != "int" )
      throw unsupported();

    current_pred.par[name] = pos;
  }

  virtual void predicateExpression(AST *tree)
  {
    current_pred.tree = tree;
  }

  virtual void endPredicate()
  {
    preds[current_pred_name] = current_pred;
    posters[current_pred_name] = &minicsp_callback::post_expression;
  }

  virtual void beginConstraint(const string & name, int idConstr,
                               int arity,
                               const string & reference,
                               CSPDefinitionType type, int id,
                               const ASTList &scope)
  {
    constraintReference=reference;
    current_scope = ast2vararray(scope);
  }

  virtual void constraintParameters(const ASTList &args)
  {
    if( posters.find(constraintReference) == posters.end() ) {
      cout << "could not find poster for reference " << constraintReference << "\n";
      throw unsupported();
    }
    (this->*posters[constraintReference])(constraintReference, args);
  }

  /********************************************************************/


  /**
   * signal the end of parsing
   */
  virtual void endInstance()
  {
  }

  void print_solution()
  {
    map<string, cspvar>::const_iterator i, b = vars.begin(), e = vars.end();
    for(i = b; i != e; ++i) {
      if( i != b ) cout << ", ";
      cout << cspvar_printer(_solver, i->second) << "=" << _solver.cspModelValue(i->second);
    }
    cout << "\n";
  }
};

} // namespace

#endif

