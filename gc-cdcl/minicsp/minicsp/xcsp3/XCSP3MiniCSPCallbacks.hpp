#ifndef _MINICSP_XCSP3_CALLBAKCS_H
#define _MINICSP_XCSP3_CALLBAKCS_H

#include <map>
#include <variant>
#include "XCSP3CoreCallbacks.h"
#include "XCSP3Variable.h"
#include "minicsp/core/cons.hpp"
#include "Tree.hpp"


namespace XCSP3Core {
    using namespace minicsp;

    struct ConstantObjective {};
    struct MinimizeObjective { cspvar var; };
    struct MaximizeObjective { cspvar var; };

    enum ObjectiveSense { SENSE_MINIMIZE, SENSE_MAXIMIZE };

    class XCSP3MiniCSPCallbacks : public XCSP3CoreCallbacks {


    public:
        Solver &solver;
        map<string, cspvar> tocspvars;
        map<int, cspvar> constants;

        vector<vector<int>> *previousTuples;
        vector<int> previousTuplesSize1;

        std::variant<ConstantObjective, MinimizeObjective, MaximizeObjective>
            objective{ConstantObjective{}};

        XCSP3MiniCSPCallbacks(Solver &s) : XCSP3CoreCallbacks(), solver(s) {
            recognizeSpecialCountCases = false;
            recognizeNValuesCases = false;
            recognizeSpecialIntensionCases = false;
        }


        void print_solution() {
            map<string, cspvar>::const_iterator i, b = tocspvars.begin(), e = tocspvars.end();
            cout << "\nv <instantiation type='solution'>\nv <list>";
            for(i = b; i != e; ++i)
                cout << i->first << " ";
            cout << "v </list>\nv <values>";
            for(i = b; i != e; ++i)
                cout << solver.cspModelValue(i->second) << " ";

            cout << "\nv </values>\n</instantiation>\n";
        }

        // ---------------------------- StartInstance -------------------------------

        void beginInstance(InstanceType type) override {}

        // ---------------------------- XVariable -> minicsp variables -------------------------------

        vector<cspvar> xvars2cspvars(const vector<XVariable *> &xvars) {
          vector<cspvar> v(xvars.size());
          std::transform(xvars.begin(), xvars.end(), v.begin(),
                         [this](XVariable *x) { return tocspvars[x->id]; });
          return v;
        }

        cspvar constant(int c) {
          if (constants.find(c) == constants.end())
            constants[c] = solver.newCSPVar(c, c);
          return constants[c];
        }

        // ---------------------------- VARIABLES ------------------------------------------

        void buildVariableInteger(string id, int minValue, int maxValue) override {
            cspvar x = solver.newCSPVar(minValue, maxValue);
            solver.setCSPVarName(x, id);
            tocspvars[id] = x;
        }


        void buildVariableInteger(string id, vector<int> &values) override {
            for(size_t i = 0; i < values.size() - 1; i++)
                if(values[i + 1] == values[i])
                    throw runtime_error("Probem : domain with identical value: " + id);

            cspvar x = solver.newCSPVar(values[0], values.back());
            int v = values[0];
            for(size_t i = 1; i < values.size(); i++) {
                v++;
                while(v != values[i])
                    x.remove(solver, v++, NO_REASON);
            }
            solver.setCSPVarName(x, id);
            tocspvars[id] = x;
        }

        // -------------------- OBJECTIVE --------------------
        void buildObjectiveMinimizeExpression(string expr) override {
            objective = MinimizeObjective{varFromExpression(expr)};
        }

        void buildObjectiveMaximizeExpression(string expr) override {
            objective = MaximizeObjective{varFromExpression(expr)};
        }

        void buildObjectiveMinimizeVariable(XVariable *x) override {
            objective = MinimizeObjective{tocspvars[x->id]};
        }

        void buildObjectiveMaximizeVariable(XVariable *x) override {
            objective = MaximizeObjective{tocspvars[x->id]};
        }

        // build up a variable that holds the value of the special
        // form objective and return that value. We get coefs by
        // value, because we may modify them
        cspvar special_form_objective(ObjectiveSense sense,
                                      ExpressionObjective type,
                                      vector<XVariable *> const &list,
                                      vector<int> coefs) {
          if (list.size() != coefs.size())
            throw runtime_error("number of coefficients does not match "
                                "number of variables in objective");
          cspvar rv;
          auto xs = xvars2cspvars(list);
          switch (type) {
          case EXPRESSION_O:
            throw runtime_error(
                "expression objective used as special form objective");
          case SUM_O: {
            int omin{0}, omax{0};
            for (size_t i = 0; i != list.size(); ++i) {
              auto x = xs[i];
              if (coefs[i] > 0) {
                omin += coefs[i] * x.min(solver);
                omax += coefs[i] * x.max(solver);
              } else if (coefs[i] < 0) {
                omin += coefs[i] * x.max(solver);
                omax += coefs[i] * x.min(solver);
              }
            }
            rv = solver.newCSPVar(omin, omax);
            xs.push_back(rv);
            coefs.push_back(-1);
            if (sense == SENSE_MAXIMIZE) {
              for (auto &c : coefs)
                c = -c;
            }
            post_lin_leq(solver, xs, coefs, 0);

          } break;
          case PRODUCT_O:
            throw runtime_error(
                "product with special form objective not supported");
          case MINIMUM_O: // fallthrough
          case MAXIMUM_O: {
            int omin{INT_MAX}, omax{INT_MIN};
            for (size_t i = 0; i != list.size(); ++i) {
              auto x = xs[i];
              if (coefs[i] > 0) {
                omin = std::min(omin, coefs[i] * x.min(solver));
                omax = std::max(omax, coefs[i] * x.max(solver));
              } else if (coefs[i] < 0) {
                omax = std::min(omin, coefs[i] * x.min(solver));
                omin = std::max(omax, coefs[i] * x.max(solver));
              }
            }
            rv = solver.newCSPVar(omin, omax);
            for (size_t i = 0; i != list.size(); ++i) {
              // need view variables here
              assert(0);
            }
          } break;
          case NVALUES_O:
            throw runtime_error(
                "nvalues special form objective not supported");
          case LEX_O:
            throw runtime_error(
                "lexicographic objective not supported");
          }
          return rv;
        }

        void buildObjectiveMinimize(ExpressionObjective type,
                                    vector<XVariable *> &list,
                                    vector<int> &coefs) override {
          cspvar obj{special_form_objective(SENSE_MINIMIZE, type, list, coefs)};
          objective = MinimizeObjective{obj};
        }

        void buildObjectiveMaximize(ExpressionObjective type,
                                    vector<XVariable *> &list,
                                    vector<int> &coefs) override {
          cspvar obj{special_form_objective(SENSE_MAXIMIZE, type, list, coefs)};
          objective = MaximizeObjective{obj};
        }

        // ---------------------------- EXTENSION ------------------------------------------

        void buildConstraintExtension(string id, vector<XVariable *> list,
                                      vector<vector<int>> &tuples, bool support,
                                      bool hasStar) override {
          if (hasStar) {
            // this should be a no-op transformation, but we map from
            // STAR to STAR_CONSTANT in case either of these changes
            // in the future
            for (auto &t : tuples) {
              for (auto &v : t)
                if (v == STAR)
                  v = STAR_CONSTANT;
            }
          }
          previousTuples = &tuples;
          if (support)
            post_positive_table(solver, xvars2cspvars(list), tuples);
          else
            post_negative_table(solver, xvars2cspvars(list), tuples);
        }


        void buildConstraintExtension(string id, XVariable *variable, vector<int> &tuple, bool support, bool hasStar) override {
            if(hasStar)
                throw runtime_error("* not supported in extensional constraints");
            previousTuplesSize1 = tuple;

            cspvar x = tocspvars[variable->id];
            if(support) {
                for(int v = x.min(solver); v <= x.max(solver); v++)
                    if(x.indomain(solver, v)) {
                        std::vector<int>::iterator it = std::find(tuple.begin(), tuple.end(), v);
                        if(it == tuple.end())
                            x.remove(solver, v, NO_REASON);
                    }
            } else {
                for(int v : tuple)
                    if(x.indomain(solver, v)) x.remove(solver, v, NO_REASON);
            }
        }


        void buildConstraintExtensionAs(string id, vector<XVariable *> list, bool support, bool hasStar) override {
            if(list.size() == 1)
                buildConstraintExtension(id, list[0], previousTuplesSize1, support, hasStar);
            else
                buildConstraintExtension(id, list, *previousTuples, support, hasStar);
        }


        // ---------------------------- INTENSION + PRIMITIVES ------------------------------------------

        cspvar postExpression(Node *n, bool isRoot = false);

        cspvar varFromExpression(string expr) {
            Tree tree(expr);
            return postExpression(tree.root, true);
        }

        void buildConstraintIntension(string id, string expr) override {
            Tree tree(expr);
            postExpression(tree.root, true);
        }

        void buildConstraintIntension(string id, Tree* tree) override {
            postExpression(tree->root, true);
        }

        // ---------------------------- LANGUAGES ------------------------------------------

        void
        buildConstraintRegular(string id, vector<XVariable *> &list, string st, vector<string> &final, vector<XTransition> &transitions) override {

            map<string, size_t> states;
            size_t current = 1;
            for(XTransition xt : transitions) {
                if(states.find(xt.from) == states.end())
                    states[xt.from] = current++;
                if(states.find(xt.to) == states.end())
                    states[xt.to] = current++;
            }

            vector<regular::transition> minitransitions;
            for(XTransition xt : transitions) {
                regular::transition t(states[xt.from], xt.val, states[xt.to]);
                minitransitions.push_back(t);
            }
            set<int> finals;
            for(string f : final)
                finals.insert(states[f]);
            regular::automaton aut(minitransitions, states[st], finals);
            post_regular(solver, xvars2cspvars(list), aut);
        }

        // ---------------------------- ALLDIFF ALLEQUAL ------------------------------------------

        void buildConstraintAlldifferent(string id, vector<XVariable *> &list) override {
            post_alldiff(solver, xvars2cspvars(list));
        }


        void buildConstraintAllEqual(string id, vector<XVariable *> &list) override {
            vector<cspvar> vars = xvars2cspvars(list);
            for(size_t i = 0; i < vars.size() - 1; i++)
                post_eq(solver, vars[i], vars[i + 1], 0);
        }

        // ---------------------------- COUNT ------------------------------------------
        void buildConstraintCount(string id, vector<XVariable *> &list,
                                  vector<int> &vals, XCondition &xc) override {
          auto xs = xvars2cspvars(list);

          vector<Var> vars;
          vector<Lit> ps;
          for (auto x : xs) {
            if (vals.size() == 1) {
              vars.push_back(x.eqi(solver, vals[0]));
              continue;
            }
            Var var = solver.newVar();
            vars.push_back(var);
            ps.clear();
            ps.push_back(~Lit(var));
            for (auto v : vals) {
              if (!x.indomain(solver, v))
                continue;
              ps.push_back(x.r_neq(solver, v));
              solver.addClause(vector<Lit>{x.r_eq(solver, v), Lit(var)});
            }
            solver.addClause(ps);
          }

          vector<int> poscoeff(vars.size(), 1), negcoeff(vars.size(), -1);
          if (xc.op == IN) {
            if (xc.operandType == INTERVAL) {
              post_pb(solver, vars, poscoeff, xc.min);
              post_pb(solver, vars, negcoeff, -xc.max);
            } else {
              auto rhs = tocspvars[xc.var];
              post_pb(solver, vars, poscoeff, 0, rhs);
              assert(0);
              post_pb(solver, vars, negcoeff, 0, rhs /* XXX: but negated */);
            }
          } else {
            int correction{0};
            if (xc.op == LT || xc.op == GT)
              correction = 1;
            if (xc.op == LT || xc.op == LE || xc.op == EQ)
              post_pb(solver, vars, negcoeff, -(xc.val + correction));
            if (xc.op == GT || xc.op == GE || xc.op == EQ)
              post_pb(solver, vars, poscoeff, xc.val + correction);
            if (xc.op == NE) {
              // build an automaton
              throw runtime_error("count() != constant not supported");
            }
          }
        }

        // ---------------------------- ORDERED ------------------------------------------

        void buildConstraintOrdered(string id, vector<XVariable *> &list, OrderType order) override {
            vector<cspvar> vars = xvars2cspvars(list);
            for(size_t i = 0; i < vars.size() - 1; i++) {
                if(order == LE)
                    post_leq(solver, vars[i], vars[i + 1], 0);
                if(order == LT)
                    post_less(solver, vars[i], vars[i + 1], 0);
                if(order == GE)
                    post_leq(solver, vars[i + 1], vars[i], 0);
                if(order == GT)
                    post_less(solver, vars[i + 1], vars[i], 0);
            }
        }


        void buildConstraintLex(string id, vector<vector<XVariable *>> &lists, OrderType order) override {
            vector<cspvar> vars1, vars2;
            for(size_t i = 0; i < lists.size() - 1; i++) {
                vars1 = xvars2cspvars(lists[i]);
                vars2 = xvars2cspvars(lists[i + 1]);
                if(order == LE)
                    post_lex_leq(solver, vars1, vars2);
                if(order == LT)
                    post_lex_less(solver, vars1, vars2);
                if(order == GE)
                    post_lex_leq(solver, vars2, vars1);
                if(order == GT)
                    post_lex_less(solver, vars2, vars1);
            }
        }


        void buildConstraintLexMatrix(string id, vector<vector<XVariable *>> &matrix, OrderType order) override {
            vector<cspvar> vars1, vars2;
            // lines
            buildConstraintLex(id, matrix, order);

            //columns
            vector<vector<XVariable *>> tmatrix;
            for(size_t i = 0 ; i < matrix[0].size() ; i++) {
                vector<XVariable *>tmp;
                for(size_t j = 0 ; j < matrix.size() ; j++)
                    tmp.push_back(matrix[j][i]);
                tmatrix.push_back(tmp);
            }
            buildConstraintLex(id, tmatrix, order);
        }


        // ---------------------------- SUM ------------------------------------------

        void postSum(string id, vector<cspvar> &list, vector<int> &coefs, XCondition &xc) {
            xc.val = -xc.val;
            switch(xc.op) {
                case EQ:
                    post_lin_eq(solver, list, coefs, xc.val);
                    break;
                case NE:
                    post_lin_neq(solver, list, coefs, xc.val);
                    break;
                case GE:
                    for(size_t i = 0; i != coefs.size(); ++i) coefs[i] = -coefs[i];
                    post_lin_leq(solver, list, coefs, -xc.val);
                    break;
                case GT:
                    for(size_t i = 0; i != coefs.size(); ++i) coefs[i] = -coefs[i];
                    post_lin_less(solver, list, coefs, -xc.val);
                    break;
                case LE:
                    post_lin_leq(solver, list, coefs, xc.val);
                    break;
                case LT:
                    post_lin_less(solver, list, coefs, xc.val);
                    break;
                default:
                    throw runtime_error("this sum is not supported");
            }
        }


        void buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coeffs, XCondition &xc) override {
            vector<cspvar> variables = xvars2cspvars(list);
            if(xc.operandType == VARIABLE) {
                xc.operandType = INTEGER;
                xc.val = 0;
                variables.push_back(tocspvars[xc.var]);
                coeffs.push_back(-1);
            }
            if(xc.op != IN) {
                postSum(id, variables, coeffs, xc);
                return;
            }
            // Intervals
            xc.op = GE;
            xc.val = xc.min;
            postSum(id, variables, coeffs, xc);
            xc.op = LE;
            xc.val = xc.max;
            postSum(id, variables, coeffs, xc);
        }


        void buildConstraintSum(string id, vector<XVariable *> &list, XCondition &xc) override {
            vector<int> coeffs;
            coeffs.assign(list.size(), 1);
            buildConstraintSum(id, list, coeffs, xc);
        }


        // Simulate scalar sum with instension constraint.

        void buildConstraintSum(string id, vector<XVariable *> &list, vector<XVariable *> &coeffs, XCondition &xc) override {
            string tmp = "add(";
            assert(list.size() == coeffs.size());
            for(size_t i = 0 ; i < list.size() ; i++) {
                tmp = tmp + "mul(" + list[i]->id + "," + coeffs[i]->id + ")";
                if(i < list.size() - 1) tmp = tmp + ",";
            }
            if(xc.operandType == VARIABLE) {
                xc.operandType = INTEGER;
                xc.val = 0;
                tmp = tmp + ",neg(" + xc.var + ")";
            }
            tmp = tmp + ")";
            if(xc.op != IN) {
                if(xc.op == EQ) tmp = "eq(" + tmp;
                if(xc.op == NE) tmp = "ne(" + tmp;
                if(xc.op == LE) tmp = "le(" + tmp;
                if(xc.op == LT) tmp = "lt(" + tmp;
                if(xc.op == GE) tmp = "ge(" + tmp;
                if(xc.op == GT) tmp = "gt(" + tmp;
                tmp = tmp + "," + std::to_string(xc.val) + ")";
                buildConstraintIntension(id,tmp);
                return;
            }

            // Intervals
            buildConstraintIntension(id,"ge("+tmp+","+std::to_string(xc.min));
            buildConstraintIntension(id,"le("+tmp+","+std::to_string(xc.max));
        }



        // ---------------------------- AtMostNValues ------------------------------------------

        void buildConstraintNValues(string id, vector<XVariable *> &list, XCondition &xc) override {
            if(xc.op != LE)
                throw runtime_error("nValues is only supported with atMost");
            vector<cspvar> vars = xvars2cspvars(list);
            cspvar n;
            switch(xc.operandType) {
                case INTEGER :
                    n = constant(xc.val);
                    break;
                case VARIABLE :
                    n = tocspvars[xc.var];
                    break;
                default :
                    throw runtime_error("nValues with interval is not yet supported");
            }
            post_atmostnvalue(solver, vars, n);
        }

        // ---------------------------- MIN/MAX ------------------------------------------
        // Use expression to do that, avoid duplication of code....

        string opToString(OrderType o) {
            if(o == LE) return "le(";
            if(o == LT) return "lt(";
            if(o == GE) return "ge(";
            if(o == GT) return "gt(";
            if(o == EQ) return "eq(";
            if(o == NE) return "ne(";
            throw runtime_error("Strange...");
        }


        void createMinMaxExpression(string minmax, vector<XVariable *> &list, XCondition &xc) {
            string tmp = minmax + "(";
            for(XVariable *xv : list)
                tmp += xv->id + (xv == list.back() ? "" : ",");
            tmp += ")";
            switch(xc.operandType) {
                case INTEGER :
                    buildConstraintIntension("", opToString(xc.op) + tmp + "," + to_string(xc.val) + ")");
                    break;
                case VARIABLE :
                    buildConstraintIntension("", opToString(xc.op) + tmp + "," + xc.var + ")");
                    break;
                case INTERVAL :
                    buildConstraintIntension("", "ge(" + tmp + "," + to_string(xc.min) + ")");
                    buildConstraintIntension("", "le(" + tmp + "," + to_string(xc.max) + ")");
                    break;
            }


        }


        virtual void buildConstraintMinimum(string id, vector<XVariable *> &list, XCondition &xc) override {
            createMinMaxExpression("min", list, xc);
        }


        virtual void buildConstraintMaximum(string id, vector<XVariable *> &list, XCondition &xc) override {
            createMinMaxExpression("max", list, xc);
        }



        // ---------------------------- ELEMENT ------------------------------------------

        void
        buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int value) override {
            if(rank != ANY)
                throw runtime_error("Basic element is only supported");
            post_element(solver, constant(value), tocspvars[index->id], xvars2cspvars(list), startIndex);
        }


        void
        buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) override {
            if(rank != ANY)
                throw runtime_error("Basic element is only supported");
            post_element(solver, tocspvars[value->id], tocspvars[index->id], xvars2cspvars(list), startIndex);
        }

        // ---------------------------- CHANNEL -------------------------------
        void buildConstraintChannel(string id, vector<XVariable *> &list,
                                    int startIndex) override {
          auto vs = xvars2cspvars(list);

          for (size_t i = 0; i != vs.size(); ++i) {
            vs[i].setmin(solver, startIndex, NO_REASON);
            vs[i].setmax(solver, vs.size() + startIndex, NO_REASON);
            for (size_t j = 0; j != vs.size(); ++j) {
              if (!vs[i].indomain(solver, j + startIndex))
                continue;
              else if (!vs[j].indomain(solver, i + startIndex))
                vs[i].remove(solver, j, NO_REASON);
              else {
                solver.addClause(
                    std::vector<Lit>({vs[i].r_eq(solver, j + startIndex),
                                      vs[j].e_eq(solver, i + startIndex)}));
              }
            }
          }
        }

        void buildConstraintChannel(string id, vector<XVariable *> &list1,
                                    int startIndex1, vector<XVariable *> &list2,
                                    int startIndex2) override {
          auto xs = xvars2cspvars(list1);
          auto ys = xvars2cspvars(list2);

          for (size_t i = 0; i != xs.size(); ++i) {
            xs[i].setmin(solver, startIndex2, NO_REASON);
            xs[i].setmax(solver, ys.size() + startIndex2, NO_REASON);
          }
          for (size_t i = 0; i != ys.size(); ++i) {
            ys[i].setmin(solver, startIndex2, NO_REASON);
            ys[i].setmax(solver, xs.size() + startIndex1, NO_REASON);
          }

          for (size_t i = 0; i != xs.size(); ++i) {
            for (size_t j = 0; j != ys.size(); ++j) {
              if (!xs[i].indomain(solver, j + startIndex2))
                continue;
              else if (!ys[j].indomain(solver, i + startIndex1))
                xs[i].remove(solver, j, NO_REASON);
              else {
                solver.addClause(
                    std::vector<Lit>({xs[i].r_eq(solver, j + startIndex2),
                                      ys[j].e_eq(solver, i + startIndex1)}));
              }
            }
          }
        }

        void buildConstraintChannel(string id, vector<XVariable *> &list,
                                    int startIndex, XVariable *value) override {
          auto bs = xvars2cspvars(list);
          auto x = tocspvars[value->id];

          x.setmin(solver, startIndex, NO_REASON);
          x.setmax(solver, startIndex + bs.size(), NO_REASON);
          for (size_t i = 0; i != bs.size(); ++i) {
            if (!x.indomain(solver, i + startIndex))
              bs[i].setmax(solver, 0, NO_REASON);
            else if (bs[i].max(solver) == 0)
              x.remove(solver, i + startIndex, NO_REASON);
            else if (bs[i].min(solver) == 1)
              x.assign(solver, i + startIndex, NO_REASON);
            else {
              solver.addClause(std::vector<Lit>{x.r_eq(solver, i + startIndex),
                                                bs[i].e_eq(solver, 1)});
              solver.addClause(std::vector<Lit>{x.r_neq(solver, i + startIndex),
                                                bs[i].e_eq(solver, 0)});
            }
          }
        }

        // ---------------------------- Instantiation -----------------------
        void buildConstraintInstantiation(string id, vector<XVariable *> &list,
                                          vector<int> &values) override {
          for (size_t i = 0; i < list.size(); i++)
            tocspvars[list[i]->id].assign(solver, values[i], NO_REASON);
        }
    };

    // -----------------------------------------------------------------------
    // ---------------------------- INTENSIONAL : POST EXPRESSION ! ----------
    // -----------------------------------------------------------------------
    cspvar XCSP3MiniCSPCallbacks::postExpression(Node *n, bool root) {
        cspvar rv;
        if(n->type == OVAR) {
            assert(!root);
            NodeVariable *nv = (NodeVariable *) n;
            return tocspvars[nv->var];
        }

        if(n->type == ODECIMAL) {
            assert(!root);
            NodeConstant *nc = (NodeConstant *) n;
            return constant(nc->val);

        }

        NodeOperator *fn = (NodeOperator *) n;

        if(fn->type == OEQ) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_eq(solver, x1, x2, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_eq_re(solver, x1, x2, 0, rv);
            }
        }

        if(fn->type == ONE) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_neq(solver, x1, x2, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_neq_re(solver, x1, x2, 0, rv);
            }
        }

        if(fn->type == OGE) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_leq(solver, x2, x1, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_geq_re(solver, x1, x2, 0, rv);
            }
        }

        if(fn->type == OGT) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_less(solver, x2, x1, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_gt_re(solver, x1, x2, 0, rv);
            }
        }


        if(fn->type == OLE) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_leq(solver, x1, x2, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_leq_re(solver, x1, x2, 0, rv);
            }
        }

        if(fn->type == OLT) {
            cspvar x1 = postExpression(fn->parameters[0]);
            cspvar x2 = postExpression(fn->parameters[1]);
            if(root)
                post_less(solver, x1, x2, 0);
            else {
                rv = solver.newCSPVar(0, 1);
                post_less_re(solver, x1, x2, 0, rv);
            }
        }

        if(fn->type == OIMP) { // IMP(X,Y) = NOT X OR Y
            NodeOperator *tmp = new NodeNot();
            tmp->addParameter(fn->parameters[0]);
            fn->parameters[0] = tmp;
            fn->type = OOR;
        }

        if(fn->type == OOR) {
            if(root) {
                vec<Lit> ps;
                for(size_t i = 0; i != fn->parameters.size(); ++i) {
                    cspvar arg = postExpression(fn->parameters[i]);
                    ps.push(arg.r_eq(solver, 0));
                }
                solver.addClause(ps);
            } else {
                vec<Lit> ps;
                rv = solver.newCSPVar(0, 1);
                ps.push(rv.e_eq(solver, 0));
                for(size_t i = 0; i != fn->parameters.size(); ++i) {
                    cspvar arg = postExpression(fn->parameters[i]);

                    vec<Lit> ps1;
                    ps1.push(rv.r_eq(solver, 0));
                    ps1.push(arg.e_eq(solver, 0));
                    solver.addClause(ps1);

                    ps.push(arg.r_eq(solver, 0));
                }
                solver.addClause(ps);
            }
        }

        if(fn->type == OAND) {
            if(root) {
                for(size_t i = 0; i != fn->parameters.size(); ++i)
                    postExpression(fn->parameters[i]);
            } else {
                vec<Lit> ps;
                rv = solver.newCSPVar(0, 1);
                ps.push(rv.r_eq(solver, 0));
                for(size_t i = 0; i != fn->parameters.size(); ++i) {
                    cspvar arg = postExpression(fn->parameters[i]);

                    vec<Lit> ps1;
                    ps1.push(rv.e_eq(solver, 0));
                    ps1.push(arg.r_eq(solver, 0));
                    solver.addClause(ps1);

                    ps.push(arg.e_eq(solver, 0));
                }
                solver.addClause(ps);
            }
        }

        if(fn->type == ONOT) {
            if(root) {
                cspvar x = postExpression(fn->parameters[0]);
                DO_OR_THROW(x.assign(solver, 0, NO_REASON));
            } else {
                vec<Lit> ps;
                rv = solver.newCSPVar(0, 1);
                cspvar arg = postExpression(fn->parameters[0]);
                vec<Lit> ps1, ps2;
                ps1.push(rv.r_eq(solver, 0));
                ps1.push(arg.e_neq(solver, 0));
                ps2.push(rv.r_neq(solver, 0));
                ps2.push(arg.e_eq(solver, 0));
                solver.addClause(ps1);
                solver.addClause(ps2);
            }
        }
        if(fn->type == OIFF) {
            if(root) {
                cspvar arg1 = postExpression(fn->parameters[0]),
                        arg2 = postExpression(fn->parameters[1]);
                vec<Lit> ps1, ps2;
                ps1.push(arg1.r_eq(solver, 0));
                ps1.push(arg2.e_eq(solver, 0));
                ps2.push(arg1.r_neq(solver, 0));
                ps2.push(arg2.e_neq(solver, 0));
                solver.addClause(ps1);
                solver.addClause(ps2);
            } else {
                vec<Lit> ps;
                rv = solver.newCSPVar(0, 1);
                cspvar arg1 = postExpression(fn->parameters[0]),
                        arg2 = postExpression(fn->parameters[1]);
                vec<Lit> ps1, ps2, ps3, ps4;
                ps1.push(rv.r_neq(solver, 0));
                ps1.push(arg1.r_eq(solver, 0));
                ps1.push(arg2.e_eq(solver, 0));
                ps2.push(rv.r_neq(solver, 0));
                ps2.push(arg1.r_neq(solver, 0));
                ps2.push(arg2.e_neq(solver, 0));

                ps3.push(rv.r_eq(solver, 0));
                ps3.push(arg1.r_eq(solver, 0));
                ps3.push(arg2.e_neq(solver, 0));
                ps4.push(rv.r_eq(solver, 0));
                ps4.push(arg1.r_neq(solver, 0));
                ps4.push(arg2.e_eq(solver, 0));
                solver.addClause(ps1);
                solver.addClause(ps2);
                solver.addClause(ps3);
                solver.addClause(ps4);
            }
        }

        if(fn->type == OXOR) {
            if(root) {
                cspvar arg1 = postExpression(fn->parameters[0]),
                        arg2 = postExpression(fn->parameters[1]);
                vec<Lit> ps1, ps2;
                ps1.push(arg1.r_eq(solver, 0));
                ps1.push(arg2.e_neq(solver, 0));
                ps2.push(arg1.r_neq(solver, 0));
                ps2.push(arg2.e_eq(solver, 0));
                solver.addClause(ps1);
                solver.addClause(ps2);
            } else {
                vec<Lit> ps;
                rv = solver.newCSPVar(0, 1);
                cspvar arg1 = postExpression(fn->parameters[0]),
                        arg2 = postExpression(fn->parameters[1]);
                vec<Lit> ps1, ps2, ps3, ps4;
                ps1.push(rv.r_neq(solver, 0));
                ps1.push(arg1.r_eq(solver, 0));
                ps1.push(arg2.e_neq(solver, 0));
                ps2.push(rv.r_neq(solver, 0));
                ps2.push(arg1.r_neq(solver, 0));
                ps2.push(arg2.e_eq(solver, 0));

                ps3.push(rv.r_eq(solver, 0));
                ps3.push(arg1.r_eq(solver, 0));
                ps3.push(arg2.e_eq(solver, 0));
                ps4.push(rv.r_eq(solver, 0));
                ps4.push(arg1.r_neq(solver, 0));
                ps4.push(arg2.e_neq(solver, 0));
                solver.addClause(ps1);
                solver.addClause(ps2);
                solver.addClause(ps3);
                solver.addClause(ps4);
            }
        }

        // function stuff
        if(fn->type == ONEG) {
            assert(!root);
            cspvar arg = postExpression(fn->parameters[0]);
            rv = solver.newCSPVar(-arg.max(solver), -arg.min(solver));
            post_neg(solver, arg, rv, 0);
        }

        if(fn->type == OABS) {
            assert(!root);
            cspvar arg = postExpression(fn->parameters[0]);
            rv = solver.newCSPVar(0, max(abs(arg.min(solver)), abs(arg.max(solver))));
            post_abs(solver, arg, rv, 0);
        }

        if(fn->type == OSUB) {
            assert(!root);
            cspvar arg1 = postExpression(fn->parameters[0]);
            cspvar arg2 = postExpression(fn->parameters[1]);
            int min = arg1.min(solver) - arg2.max(solver);
            int max = arg1.max(solver) - arg2.min(solver);
            rv = solver.newCSPVar(min, max);
            vector<int> w(3);
            vector<cspvar> v(3);
            w[0] = 1;
            v[0] = rv;
            w[1] = -1;
            v[1] = arg1;
            w[2] = 1;
            v[2] = arg2;
            post_lin_eq(solver, v, w, 0);
        }

        if(fn->type == ODIST) { //Simulate DIST(X,Y) = ABS(SUB(X,Y))
            fn->type = OSUB;
            // Copy of ABS Op (Above)
            cspvar arg = postExpression(fn); // call on same node
            rv = solver.newCSPVar(0, max(abs(arg.min(solver)), abs(arg.max(solver))));
            post_abs(solver, arg, rv, 0);
            fn->type = ODIST; //Be careful, depend the order,  but OSUB can be called two times otherwise :(
        }

        if(fn->type == OADD) {
            assert(!root);
            int min = 0, max = 0;
            vector<int> w;
            vector<cspvar> v;
            for(size_t q = 0; q != fn->parameters.size(); ++q) {
                cspvar arg = postExpression(fn->parameters[q]);
                w.push_back(-1);
                v.push_back(arg);
                min += arg.min(solver);
                max += arg.max(solver);
            }
            rv = solver.newCSPVar(min, max);
            v.push_back(rv);
            w.push_back(1);
            post_lin_eq(solver, v, w, 0);
        }
        if(fn->type == OMUL) {
            assert(!root);
            cspvar arg0 = postExpression(fn->parameters[0]);
            for(size_t q = 1; q != fn->parameters.size(); ++q) {
                cspvar arg1 = postExpression(fn->parameters[q]);
                int minv = min(min(arg0.min(solver) * arg1.min(solver),
                                   arg0.min(solver) * arg1.max(solver)),
                               min(arg0.max(solver) * arg1.min(solver),
                                   arg0.max(solver) * arg1.max(solver))),
                        maxv = max(max(arg0.min(solver) * arg1.min(solver),
                                       arg0.min(solver) * arg1.max(solver)),
                                   max(arg0.max(solver) * arg1.min(solver),
                                       arg0.max(solver) * arg1.max(solver)));
                cspvar res = solver.newCSPVar(minv, maxv);
                post_mult(solver, res, arg0, arg1);
                arg0 = res;
            }
            rv = arg0;
        }

        if(fn->type == OMIN) {
            assert(!root);
            cspvar arg0 = postExpression(fn->parameters[0]);
            for(size_t q = 1; q != fn->parameters.size(); ++q) {
                cspvar arg1 = postExpression(fn->parameters[q]);
                int minv = min(arg0.min(solver), arg1.min(solver)),
                        maxv = min(arg0.max(solver), arg1.max(solver));
                cspvar res = solver.newCSPVar(minv, maxv);
                post_min(solver, res, arg0, arg1);
                arg0 = res;
            }
            rv = arg0;
        }
        if(fn->type == OMAX) {
            assert(!root);
            cspvar arg0 = postExpression(fn->parameters[0]);
            for(size_t q = 1; q != fn->parameters.size(); ++q) {
                cspvar arg1 = postExpression(fn->parameters[q]);
                int minv = max(arg0.min(solver), arg1.min(solver)),
                        maxv = max(arg0.max(solver), arg1.max(solver));
                cspvar res = solver.newCSPVar(minv, maxv);
                post_max(solver, res, arg0, arg1);
                arg0 = res;
            }
            rv = arg0;
        }

        if(fn->type == OIF) {
            assert(!root);
            assert(fn->parameters.size() == 3);
            cspvar argif = postExpression(fn->parameters[0]);
            cspvar arg1 = postExpression(fn->parameters[1]);
            cspvar arg2 = postExpression(fn->parameters[2]);
            int minv = min(arg1.min(solver), arg2.min(solver)),
                    maxv = max(arg1.max(solver), arg2.max(solver));
            rv = solver.newCSPVar(minv, maxv);
            post_eq_re(solver, rv, arg1, 0, argif.e_eq(solver, 1));
            post_eq_re(solver, rv, arg2, 0, argif.e_neq(solver, 1));
        }

        return rv;
    }
}
#endif //COSOCO_XCSP3MiniCSPCallbacks_H
