/****************************************************************************************[Solver.C]
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
**************************************************************************************************/

#include "Solver.h"
#include "Sort.h"
#include <cmath>
#include <algorithm>


//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :
    // Parameters: (formerly in 'SearchParams')
    var_decay(1 / 0.95), clause_decay(1 / 0.999), random_var_freq(0.02)
  , restart_first(100), restart_inc(1.5), learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // More parameters:
    //
  , expensive_ccmin  (true)
  , polarity_mode    (polarity_true)
  , verbosity        (0)

    // Statistics: (formerly in 'SolverStats')
    //
  , starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
  , clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , ok               (true)
  , cla_inc          (1)
  , var_inc          (1)
  , qhead            (0)
  , simpDB_assigns   (-1)
  , simpDB_props     (0)
  , vertexOrderHeap  (VertexOrderLt(vertexActivity))
  , order_heap       (VarOrderLt(activity))
  , random_seed      (91648253)
  , progress_estimate(0)
  , remove_satisfied (true)
  , initialized(false)
{}


Solver::~Solver()
{
    for (int i = 0; i < clauses.size(); i++) free(clauses[i]);
    for (int i = 0; i < learnts.size(); i++) free(learnts[i]);
}

struct levelOrder_lt {
  const vec<int>& level;
  bool operator () (Lit x, Lit y) {
    if (level[var(x)] == level[var(y)]) {
      if (sign(x) && !sign(y)) {
        return true;
      }
    }
    return level[var(x)] > level[var(y)];
  }
  levelOrder_lt(const vec<int>&  lvl) : level(lvl) { }
};

int Solver::transformConflict(vec<Lit>& learnt_clause,
                              int& out_btlevel,
                              Vertex& freeVertex,
                              bool& seenSigned,
                              unsigned int& colorInfo) {
  vec<vec<Vertex> > sortedVertices;
  vec<Lit> alreadyConnected;
  const Color c0 = color(learnt_clause[0]);
  for (int i = 0; i < numberOfColors; i++) {
    for (int j = (i + 1); j < numberOfColors; j++) {
      alreadyConnected.push(lit_Undef);
    }
    sortedVertices.push();
  }

  vec<bool> canSkip; //
  for (int i = 0; i < learnt_clause.size(); i++) {
    canSkip.push(false);
  }

  freeVertex = vertex(learnt_clause[0]);
  for (int i =0; i < learnt_clause.size(); i++) {
    assert(isVertex(learnt_clause[i]) && sign(learnt_clause[i]));
    for (int j = i + 1; j < learnt_clause.size(); j++) {
      Vertex v = vertex(learnt_clause[i]);
      Vertex w = vertex(learnt_clause[j]);
      Var e = mapVerticesOnMergeVariable(v,w);
      Color c1 = color(learnt_clause[i]), c2 = color(learnt_clause[j]);
      if (level[e] == 0 && value(e) == l_False && (c1 != c2)) {
#ifndef SKIPCONNECTED
        alreadyConnected[mapColors(c1, c2)] = Lit(e, false);
#endif
      } else if (level[e] == 0 && value(e) == l_True && (c1 == c2)) {
        canSkip[j] = true;
      }
    }
    if (!canSkip[i]) {
      sortedVertices[color(learnt_clause[i])].push(vertex(learnt_clause[i]));
    }
  }
  learnt_clause.clear();
  for (int i =0; i < numberOfColors; i++) {
    if (sortedVertices[i].size() == 0) {
      continue;
    }
    for (int j = 1; j < sortedVertices[i].size(); j++) {
      Vertex v = sortedVertices[i][j];
      Vertex w = sortedVertices[i][0];
      Var e = mapVerticesOnMergeVariable(v,w);
      learnt_clause.push(Lit(e, true));
      level[e] = std::max(level[mapVerticesOnColorVariable(v,i)],
                          level[mapVerticesOnColorVariable(w,i)]);
    }
    for (int j = i + 1; j < numberOfColors; j++) {
      if (sortedVertices[j] == 0
          || alreadyConnected[mapColors(i,j)] != lit_Undef) {
        continue;
      }
      Vertex v = sortedVertices[i][0];
      Vertex w = sortedVertices[j][0];
      Var e = mapVerticesOnMergeVariable(v,w);

      assert( ((level[mapVerticesOnColorVariable(v,i)] != -1
                &&level[mapVerticesOnColorVariable(w,i)] != -1))
              || ((level[mapVerticesOnColorVariable(v,j)] != -1
                  &&level[mapVerticesOnColorVariable(w,j)] != -1)));

      int l = -1;
      if (level[mapVerticesOnColorVariable(v,i)] == -1
          || level[mapVerticesOnColorVariable(w,i)] == -1) {
        l = std::max(level[mapVerticesOnColorVariable(v,j)],
                     level[mapVerticesOnColorVariable(w,j)]);
      } else if (level[mapVerticesOnColorVariable(v,j)] == -1
          || level[mapVerticesOnColorVariable(w,j)] == -1) {
        l = std::max(level[mapVerticesOnColorVariable(v,i)],
                     level[mapVerticesOnColorVariable(w,i)]);
      } else {
        l = std::min(std::max(level[mapVerticesOnColorVariable(v,i)],
                              level[mapVerticesOnColorVariable(w,i)]),
                     std::max(level[mapVerticesOnColorVariable(v,j)],
                              level[mapVerticesOnColorVariable(w,j)]));
      }
      level[e] = l;
      learnt_clause.push(Lit(e, false));
    }
  }
  sort(learnt_clause, levelOrder_lt(level));
  seenSigned  =  sortedVertices[c0].size() > 1;
  colorInfo = 0;
  if (seenSigned) {
    colorInfo   = c0;
  } else {
    unsigned int mask = 1;
    for (int i =0; i < numberOfColors; i++) {
      if (i != c0 && (sortedVertices[i].size() != 0)) {
        colorInfo = colorInfo | mask;
      }
      mask = mask << 1;
    }
  }

  if (learnt_clause.size() == 0) {
    return vertex_unsat;
  } else if (level[var(learnt_clause[0])] <= out_btlevel) {
    return vertex_reitarate;
  }
  return vertex_sat;
}

void Solver::enqueueSafeMergeClausePropagees(Clause* clause,
                                             Vertex propagatingVertex,
                                             Color c) {
  assert(clause->size() > 0);
  Lit p = (*clause)[0];
  Vertex v, w;
  factorMerge(p, v, w);
//  assert(!isAssigned(v) || !isAssigned(w));
//  assert(isAssigned(v) || isAssigned(w));
  Lit q = Lit(mapVerticesOnColorVariable(propagatingVertex, c), sign(p));
  if (value(q) == l_Undef) uncheckedEnqueue(q, clause);
}

void Solver::enqueueUnsafeMergeClausePropagees(Clause * clause,
                                               Vertex free,
                                               bool seenSigned,
                                               unsigned int colorInfo) {
  const bool b = true;
  if (seenSigned) {
    Lit p = Lit(mapVerticesOnColorVariable(free, colorInfo), b);
    uncheckedEnqueue(p, clause);
  } else {
    unsigned int mask = 1;
    colorInfo = ~colorInfo;
    for (int i =0; i < numberOfColors; i++) {
      if (colorInfo & mask) {
        Lit p = Lit(mapVerticesOnColorVariable(free, i), b);
        if (value(p) == l_Undef) {
          uncheckedEnqueue(p, clause);
        }
      }
      mask = mask << 1;
    }
  }
}


void Solver::attachMergeClause(Clause* clause) {
  Clause& c = *clause;
  assert(c.learnt() && c.size()  > 0);
  Vertex v,w;
  for (int i =0; (i < 2 && i < c.size()); i++) {
    Var q = var(c[i]);
    if (!isPresent(q)) {
      factorMerge(q, v, w);
      primaryPointers[v].push(var(c[i]));
      primaryPointers[w].push(var(c[i]));
      locations[q].first = &(primaryPointers[v].last());
      locations[q].second = &(primaryPointers[w].last());
    }
  }
  attachClause(c);
}

void Solver::detachMergeClause(Clause* clause) {
  Clause& c = *clause;
  /* Check..
   *
  Vertex v,w, v2, w2;
  Var q;
  for (int i =0; (i < 2 && i < c.size()); i++) {
    if (watches[toInt(~c[i])].size()
        + watches[toInt(c[i])].size() == 1) {
      factorMerge(c[i], v, w);

      q = primaryPointers[v].last();
      factorMerge(q, v2, w2);
      *(locations[var(c[i])].first) = q;
      primaryPointers[v].pop();
      if (v2 == v) {
        locations[q].first = locations[var(c[i])].first;
      } else {
        assert(w2 == v);
        locations[q].second = locations[var(c[i])].first;
      }

      q = primaryPointers[w].last();
      factorMerge(q, v2, w2);
      *(locations[var(c[i])].second) = q;
      primaryPointers[w].pop();
      if (v2 == w) {
        locations[q].first = locations[var(c[i])].second;
      } else {
        assert(w2 == w);
        locations[q].second = locations[var(c[i])].second;
      }

      locations[var(c[i])].first = NULL;
      locations[var(c[i])].second = NULL;
    }
  }*/
  detachClause(c);
}

//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches   .push();          // (list for positive literal)
    watches   .push();          // (list for negative literal)
    reason    .push(NULL);
    assigns   .push(toInt(l_Undef));
    level     .push(-1);
    vertexActivity  .push(0);
    activity  .push(0),
    seen      .push(0);

    polarity    .push((char)sign);
    decision_var.push(isVertex(v));

    if (isVertex(v)) {
      insertVarOrder(v);
    }
    return v;
}

void Solver::initVertexVars(int nv, int nc) {
  if (nc < 0 || nc > 32) {
    reportf("This version only supports instances up to 32 colors!!\nSorry :( exiting..\n");
    exit(0);
  }
  numberOfVertices = nv;
  numberOfColors = nc;
  initialized = true;
  for (int i =0; i < numberOfVertices; i++) {
    colorAssignedMask.push(0);
    colorNotAssignedMask.push(0);
    colorAssignement.push(-1);
    primaryPointers.push();
    primaryPointers.last().growTo(numberOfVertices + 1); // MakeSure Size is not re..
    primaryPointers.last().clear();
    vertexSeen.push(false);
  }

  for (int i =0; i < LASTMERGEVAR + 2; i++) {
    stackPos.push(LASTMERGEVAR + 1);
    locations.push_back(std::make_pair<Var*, Var*>(NULL, NULL));
  }
}

bool Solver::isInitialized() {
  return initialized;
}

bool Solver::addClause(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);

    if (!ok)
        return false;
    else{
        // Check if clause is satisfied and remove false/duplicate literals:
        sort(ps);
        Lit p; int i, j;
        for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
            if (value(ps[i]) == l_True || ps[i] == ~p)
                return true;
            else if (1 || value(ps[i]) != l_False && ps[i] != p)
                ps[j++] = p = ps[i];
        ps.shrink(i - j);
    }

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        assert(value(ps[0]) == l_Undef);
        uncheckedEnqueue(ps[0]);
        return ok = (propagate() == NULL);
    }else{
        Clause* c = Clause_new(ps, false);
        clauses.push(c);
        attachClause(*c);
    }

    return true;
}

void Solver::attachClause(Clause& c) {
    assert(c.size() > 1);
    watches[toInt(~c[0])].push(&c);
    watches[toInt(~c[1])].push(&c);
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(Clause& c) {
    assert(c.size() > 1);
    assert(find(watches[toInt(~c[0])], &c));
    assert(find(watches[toInt(~c[1])], &c));
    remove(watches[toInt(~c[0])], &c);
    remove(watches[toInt(~c[1])], &c);
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(Clause& c) {
  if (c.learnt()) {
    detachMergeClause(&c);
  } else {
    detachClause(c);
  }
  free(&c);
}


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }


// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int backTrackLevel) {
    if (decisionLevel() > backTrackLevel){
        for (int c = trail.size()-1; c >= trail_lim[backTrackLevel]; c--) {
          Var     x  = var(trail[c]);
          assigns[x] = toInt(l_Undef);
          level[x] = -1;
          stackPos[x] = nVars();
          if (isVertex(x)) {
            if (!sign(trail[c])) {
              colorAssignedMask[vertex(x)] = 0; // There can be only one color!
              colorAssignement[vertex(x)] = -1;
            } else {
              int mask = 1;
              mask = mask << color(x);
              mask = mask ^ 0xFFFFFFFF;
              colorNotAssignedMask[vertex(x)] = colorNotAssignedMask[vertex(x)] & mask;
            }
            insertVarOrder(x);
          }
        }
        qhead = trail_lim[backTrackLevel];
        trail.shrink(trail.size() - trail_lim[backTrackLevel]);
        trail_lim.shrink(trail_lim.size() - backTrackLevel);
    }
}


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit(int polarity_mode, double random_var_freq)
{
    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (toLbool(assigns[next]) == l_Undef && decision_var[next])
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || toLbool(assigns[next]) != l_Undef || !decision_var[next]) {
        if (order_heap.empty()){
            next = var_Undef;
            break;
        } else {
            next = order_heap.removeMin();
        }
    }

    bool sign = false;
    switch (polarity_mode){
    case polarity_true:  sign = false; break;
    case polarity_false: sign = true;  break;
    case polarity_user:  sign = polarity[next]; break;
    case polarity_rnd:   sign = irand(random_seed, 2); break;
    default: assert(false); }

    return next == var_Undef ? lit_Undef : Lit(next, sign);
}

bool Solver::checkPropagationLevel(Lit l, int targetLevel) {
#ifndef CHECKPROPLEVEL
  return true;
#endif
  Vertex v, w;
  factorMerge(l, v, w);
  if(!isAssigned(v) && !isAssigned(w)) {
    return false;
  }
  if (!sign(l)) {
    if(!isAssigned(v)) {
      Color c = getAssignment(w);
      Var x = mapVerticesOnColorVariable(v, c);
      Var y = mapVerticesOnColorVariable(w, c);
      return ((level[x] == targetLevel) || (level[y] == targetLevel));
    } else if (!isAssigned(w)) {
      Color c = getAssignment(v);
      Var x = mapVerticesOnColorVariable(v, c);
      Var y = mapVerticesOnColorVariable(w, c);
      return ((level[x] == targetLevel) || (level[y] == targetLevel));
    } else {
      Color c1 = getAssignment(v);
      Var x = mapVerticesOnColorVariable(v, c1);
      Var y = mapVerticesOnColorVariable(w, c1);
      bool optionA = ((level[x] == targetLevel) || (level[y] == targetLevel));

      Color c2 = getAssignment(w);
      x = mapVerticesOnColorVariable(v, c2);
      y = mapVerticesOnColorVariable(w, c2);
      bool optionB = ((level[x] == targetLevel) || (level[y] == targetLevel));
      return (optionA || optionB);
    }
  } else {
    Color c = getAssignment(w);
    Var x = mapVerticesOnColorVariable(v, c);
    Var y = mapVerticesOnColorVariable(w, c);
    return ((level[x] == targetLevel) || (level[y] == targetLevel));
  }
}

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|
|  Description:
|    Analyze conflict and produce a reason clause.
|
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|
|  Effect:
|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
|________________________________________________________________________________________________@*/

void Solver::analyze(Clause* confl, vec<Lit>& out_learnt, int& out_btlevel)
{
#ifdef SHOWTRAIL
    vec<bool> mySeen;
    for (int i =0; i < nVars(); i++) {
      mySeen.push(false);
    }
    Clause * orgConf = confl;
#endif
    if (confl->learnt()) {
      assert(checkPropagationLevel((*confl)[0], decisionLevel()));
    }


    int pathC = 0;
    Lit p     = lit_Undef;

    vec<Lit> tmp_learnt;

    // Generate conflict clause:
    //
    tmp_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;
    out_btlevel = 0;
#ifdef PRINTDEBUG
    reportf("index:%d..\n", index);
#endif
    do{
#ifdef SHOWTRAIL
      if (confl == 0) {
          dumpTrail(orgConf, mySeen);
          reportf("pathC:%d..\n", pathC);
        }
#endif
        assert(confl != NULL);          // (otherwise should be UIP)
        Clause& c = *confl;

        if (c.learnt()) {
            claBumpActivity(c);
        }

        for (int j = (p == lit_Undef || c.learnt()) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
#ifdef SHOWTRAIL
            mySeen[var(q)] = true;
#endif
            vec<Lit> myQs;
            Vertex v, w;
            Color c;
            Vertex vp = p != lit_Undef ? vertex(p) : -1;
            if (isMerge(q)) {
              factorMerge(q, v, w);
              assert((valueOfMerge(q) == l_False) || (v == vp) || (w == vp));

              if (v == vp) {
                assert(isAssigned(w));
                c = getAssignment(w);
                myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
              } else if (w== vp) {
                assert(isAssigned(v));
                c = getAssignment(v);
                myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
              } else if (!sign(q)) {
                assert(isAssigned(v)|| isAssigned(w));
                assert(colorAssignedMask[v] != colorAssignedMask[w]);
                if (!isAssigned(v)) {
                  c = getAssignment(w);
                  myQs.push(Lit(mapVerticesOnColorVariable(v,c), false));
                  assert(value(myQs.last()) == l_False);
                  myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
                  assert(value(myQs.last()) == l_False);
                } else if (!isAssigned(w)) {
                  c = getAssignment(v);
                  myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
                  assert(value(myQs.last()) == l_False);
                  myQs.push(Lit(mapVerticesOnColorVariable(w,c), false));
                  assert(value(myQs.last()) == l_False);
                } else {
                  int posV = 0, posW = 0;
                  Color cv = getAssignment(v);
                  Color cw = getAssignment(w);
                  posV = stackPos[mapVerticesOnColorVariable(v, cv)];
                  posW = stackPos[mapVerticesOnColorVariable(w, cv)];
                  if ((posV <= index) &&  (posW <= index)) {
                    myQs.push(Lit(mapVerticesOnColorVariable(v,cv), true));
                    assert(value(myQs.last()) == l_False);
                    myQs.push(Lit(mapVerticesOnColorVariable(w,cv), false));
                    assert(value(myQs.last()) == l_False);
                  } else {
                    posV = stackPos[mapVerticesOnColorVariable(v, cw)];
                    posW = stackPos[mapVerticesOnColorVariable(w, cw)];
                    if ((posV <= index) &&  (posW <= index)) {
                      myQs.push(Lit(mapVerticesOnColorVariable(v,cw), false));
                      assert(value(myQs.last()) == l_False);
                      myQs.push(Lit(mapVerticesOnColorVariable(w,cw), true));
                      assert(value(myQs.last()) == l_False);
                    } else {
                      posV = stackPos[mapVerticesOnColorVariable(v, cv)];
                      posW = stackPos[mapVerticesOnColorVariable(w, cw)];
                      assert((posV <= index) &&  (posW <= index));
                      myQs.push(Lit(mapVerticesOnColorVariable(v,cv), true));
                      assert(value(myQs.last()) == l_False);
                      myQs.push(Lit(mapVerticesOnColorVariable(w,cw), true));
                      assert(value(myQs.last()) == l_False);
                    }
                  }
                }
              } else {
                assert(isAssigned(v) && isAssigned(w));
                assert((vp == v)
                       || (vp == w)
                       || (getAssignment(v) == getAssignment(w)));
                c = ( (vp != v) ? getAssignment(v) : getAssignment(w));
                if (v != vp) {
                  myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
                  assert(value(myQs.last()) == l_False);
                }
                if (w != vp) {
                  myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
                  assert(value(myQs.last()) == l_False);
                }
              }
            } else {
              myQs.push(q);
              assert(value(myQs.last()) == l_False);
            }

            for (int i = 0; i < myQs.size(); i++) {
              q = myQs[i];
              if ( (!seen[var(q)]) && (level[var(q)] >= 0)){
#ifdef SHOWTRAIL
                if (level[var(q)] >= decisionLevel()) {
                }
                mySeen[var(q)] = true;
#endif
                vertexBumpActivity(var(q));
                seen[var(q)] = 1;
                if (level[var(q)] >= decisionLevel())
                    pathC++;
                else{
                  tmp_learnt.push(q);
                }
            }
          }
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])]) {
        }
        p     = trail[index+1];
        confl = reason[var(p)];
        seen[var(p)] = 0;
        pathC--;

    }while ((pathC > 0) || (!isVertex(p)) || (sign(p)));
    tmp_learnt[0] = ~p;

#ifdef SHOWTRAIL
    dumpTrail(orgConf, mySeen);
#endif

    // Simplify conflict clause:
    //
    int i, j;
    if (expensive_ccmin){
        uint32_t abstract_level = 0;
        for (i = 1; i < tmp_learnt.size(); i++)
            abstract_level |= abstractLevel(var(tmp_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        // TODO(Bas) shortening does not work for now!
        tmp_learnt.copyTo(analyze_toclear);
        for (i = j = 1; i < tmp_learnt.size(); i++) {
            if ( 1
                || (reason[var(tmp_learnt[i])] == NULL)
                || (!litRedundant(tmp_learnt[i], abstract_level))) {
              tmp_learnt[j++] = tmp_learnt[i];
            }
        }
    }else{
      tmp_learnt.copyTo(analyze_toclear);
        for (i = j = 1; i < tmp_learnt.size(); i++){
            Clause& c = *reason[var(tmp_learnt[i])];
            for (int k = 1; k < c.size(); k++)
                if (!seen[var(c[k])] && level[var(c[k])] > 0){
                  tmp_learnt[j++] = tmp_learnt[i];
                    break; }
        }
    }
    max_literals += tmp_learnt.size();
    tmp_learnt.shrink(i - j);
    tot_literals += tmp_learnt.size();

#ifdef SHOWORGLEARNT
    reportf("\nOrgLearnt:\n");
    printClauseFancy(tmp_learnt);
    reportf("\n\n");
#endif

    assert(level[var(p)] == decisionLevel());

    for (int i =0; i < tmp_learnt.size(); i++) {
      if (isVertex(tmp_learnt[i]) && sign(tmp_learnt[i])) {
        out_learnt.push(tmp_learnt[i]) ;
      } else {
        vec<Lit> vertexReason;
        deduceVertexAndIntroducedReason(tmp_learnt[i], vertexReason);
        for (int j =0; j < vertexReason.size(); j++) {
          out_learnt.push(vertexReason[j]);
        }
      }
    }

    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)

    sort(out_learnt, levelOrder_lt(level));

#ifdef SHOWEXTENDEGLEARNT
    reportf("\n Learnt:\n");
    printClauseFancy(out_learnt);
    reportf("\n\n");
#endif

    if (out_learnt.size() == 1) {
      out_btlevel = 0;
      reportf("Almost proven UNSAT!!!\n");
    } else {
      out_btlevel = level[var(out_learnt[1])];
    }

#ifdef SHOWTRAIL
//    dumpTrail(orgConf, mySeen);

    reportf("Learnt:\n\n");
    printClauseFancy(out_learnt);
    reportf("\n\n");
#endif

#ifdef PRINTDEBUG
    reportf("\n\nLearnt bt:%d..\n\n", out_btlevel);
    printClauseFancy(out_learnt);
    reportf("\n\n");
#endif
}


void Solver::deduceVertexAndIntroducedReason(Lit x, vec<Lit>& myReason) { // NOLINT
  assert(reason[var(x)] != NULL);
  Clause& myClause = *reason[var(x)];
  for (int i =  (myClause.learnt() ? 0 : 1) ; i < myClause.size(); i++) {
    Lit q =  (*reason[var(x)])[i];
    vec<Lit> myQs;
    Vertex v, w;
    Color c;
    if (isMerge(q)) {
      assert((reason[var(x)])->learnt());
      factorMerge(q, v, w);
      Vertex vp = vertex(x);
      if (v == vp) {
        assert(isAssigned(w));
        c = getAssignment(w);
        myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
      } else if (w == vp) {
        assert(isAssigned(v));
        c = getAssignment(v);
        myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
      } else if (!sign(q)) {
        assert(isAssigned(v) || isAssigned(w));
        assert(colorAssignedMask[v] != colorAssignedMask[w]);
        if (!isAssigned(v)) {
          c = getAssignment(w);
          myQs.push(Lit(mapVerticesOnColorVariable(v,c), false));
          assert(value(myQs.last()) == l_False);
          myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
          assert(value(myQs.last()) == l_False);
        } else if (!isAssigned(w)) {
          c = getAssignment(v);
          myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
          assert(value(myQs.last()) == l_False);
          myQs.push(Lit(mapVerticesOnColorVariable(w,c), false));
          assert(value(myQs.last()) == l_False);
        } else {
          Color cv = getAssignment(v);
          Color cw = getAssignment(w);
          int maxL = decisionLevel() + 1;

          int la;
          if (level[mapVerticesOnColorVariable(v,cv)] == -1
              || level[mapVerticesOnColorVariable(w,cw)] == -1) {
            la = maxL;
          } else {
            la  = std::max(level[mapVerticesOnColorVariable(v,cv)],
                            level[mapVerticesOnColorVariable(w,cw)]);
          }

          int lv;
          if (level[mapVerticesOnColorVariable(v,cv)] == -1
              || level[mapVerticesOnColorVariable(w,cv)] == -1) {
            lv = maxL;
          }  else {
            lv = std::max(level[mapVerticesOnColorVariable(v,cv)],
                            level[mapVerticesOnColorVariable(w,cv)]);
          }

          int lw;
          if (level[mapVerticesOnColorVariable(v,cw)] == -1
              || level[mapVerticesOnColorVariable(w,cw)] == -1) {
            lw = maxL;
          } else {
            lw = std::max(level[mapVerticesOnColorVariable(v,cw)],
                          level[mapVerticesOnColorVariable(w,cw)]);
          }

          if (la <= lv && la <= lw) {
            myQs.push(Lit(mapVerticesOnColorVariable(v,cv), true));
            assert(value(myQs.last()) == l_False);
            myQs.push(Lit(mapVerticesOnColorVariable(w,cw), true));
            assert(value(myQs.last()) == l_False);
          }
          else if (lv <= la && lv <= lw) {
            myQs.push(Lit(mapVerticesOnColorVariable(v,cv), true));
            assert(value(myQs.last()) == l_False);
            myQs.push(Lit(mapVerticesOnColorVariable(w,cv), false));
            assert(value(myQs.last()) == l_False);
          }
          else {
            myQs.push(Lit(mapVerticesOnColorVariable(v,cw), false));
            assert(value(myQs.last()) == l_False);
            myQs.push(Lit(mapVerticesOnColorVariable(w,cw), true));
            assert(value(myQs.last()) == l_False);
          }
        }
      } else {
        assert(isAssigned(v) && isAssigned(w));
        assert(getAssignment(v) == getAssignment(w));
        c =  getAssignment(v);
        if (v != vp) {
          myQs.push(Lit(mapVerticesOnColorVariable(v,c), true));
          assert(value(myQs.last()) == l_False);
        }
        if (w != vp) {
          myQs.push(Lit(mapVerticesOnColorVariable(w,c), true));
          assert(value(myQs.last()) == l_False);
        }
      }
    } else {
      myQs.push(q);
    }

    for (int i = 0; i < myQs.size(); i++) {
      if (seen[var(myQs[i])] == 1) {
        continue;
      } else {
        seen[var(myQs[i])] = 1;
        analyze_toclear.push(myQs[i]);
      }

      if ((isVertex(myQs[i]) && sign(myQs[i]))) {
        myReason.push(myQs[i]);
      } else {
        deduceVertexAndIntroducedReason(myQs[i], myReason);
      }
    }
  }
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason[var(analyze_stack.last())] != NULL);
        Clause& c = *reason[var(analyze_stack.last())]; analyze_stack.pop();

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level[var(p)] > 0){
                if (reason[var(p)] != NULL && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason[x] == NULL){
                assert(level[x] > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = *reason[x];
                for (int j = 1; j < c.size(); j++)
                    if (level[var(c[j])] > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }
    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, Clause* from)
{
    assert(value(p) == l_Undef);
    assigns [var(p)] = toInt(lbool(!sign(p)));  // <<== abstract but not uttermost effecient
    level   [var(p)] = decisionLevel();
    reason  [var(p)] = from;
    stackPos [var(p)] = trail.size();
    trail.push(p);
    if (!isVertex(p)) {
      return;
    }
    int x = 1;
    x = x << color(p);
    if (sign(p)) {
      colorNotAssignedMask[vertex(p)] = colorNotAssignedMask[vertex(p)] | x;
      assert(colorNotAssignedMask[vertex(p)] != 0);
    } else {
      colorAssignement[vertex(p)] = color(p);
      colorAssignedMask[vertex(p)] =  x;
      assert(colorAssignedMask[vertex(p)] != 0);
    };
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise NULL.
|
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
Clause* Solver::propagate() {
  Clause* confl     = NULL;
  int     num_props = 0;
  Vertex  v, w, v_, w_;

  Clause **i, **j, **end;
  Var *layerI, *layerEnd;
  vec<Vertex> vertexStack;
  do {
    while (qhead < trail.size()) {
      num_props++;
      Lit p   = trail[qhead++];  // 'p' is enqueued fact to propagate.
      if (isVertex(p) && !vertexSeen[vertex(p)]) {
        vertexSeen[vertex(p)] = true;
        vertexStack.push(vertex(p));
      }
      vec<Clause*>& ws = watches[toInt(p)];
      for (i = j = static_cast<Clause**>(ws), end = i + ws.size();  i != end;) {
        Clause& c = **i++;
        // Make sure the false literal is data[1]:
        Lit falselit = ~p;
        if (c[0] == falselit) {
          c[0] = c[1];
          c[1] = falselit;
        }
        assert(c[1] == falselit);
        // If 0th watch is true, then clause is already satisfied.
        Lit first = c[0];
        if (value(first) == l_True) {
          *j++ = &c;
        } else {
          // Look for new watch:
          for (int k = 2; k < c.size(); k++) {
            if (value(c[k]) != l_False) {
              c[1] = c[k];
              c[k] = falselit;
              watches[toInt(~c[1])].push(&c);
              goto FoundWatch;
            }
          }

          // Did not find watch -- clause is unit under assignment:
          *j++ = &c;
          if (value(first) == l_False) {
            confl = &c;
            qhead = trail.size();
            // Copy the remaining watches:
            while (i < end) {
              *j++ = *i++;
            }
          } else {
            uncheckedEnqueue(first, &c);
          }
        }
      FoundWatch: {}
      }
      ws.shrink(i - j);
    }

    for (int x =0; x < vertexStack.size(); x++) {
      if (confl != NULL) {
        continue;
      }
      const Vertex  updated = vertexStack[x];
      vec<Var> vp;
      for (int x = 0; x < primaryPointers[updated].size(); x++) {
        Var q_ = primaryPointers[updated][x];
        vp.push(q_);
        factorMerge(q_, v_, w_);
        if (v_ == updated)
          locations[q_].first = NULL;
        else
          locations[q_].second = NULL;
      }
      primaryPointers[updated].clear();
      for (layerI = static_cast<Var*>(vp), layerEnd = layerI + vp.size();
           layerI != layerEnd;
           /* empty */) {
        Var q = *layerI++;
        lbool res = valueOfMerge(q);
        if (res == l_Undef || confl != NULL) {
          int x = 0;
          factorMerge(q, v, w);
          if (locations[q].first == NULL) {
            primaryPointers[v].push(q);
            locations[q].first = &(primaryPointers[v].last());
            x++;
          }

          if (locations[q].second == NULL) {
            primaryPointers[w].push(q);
            locations[q].second = &(primaryPointers[w].last());
            x++;
          }
          if ((value(q) == l_False && isAssigned(updated) && confl == NULL)
              /*|| (value(q) == l_True)*/ ) {
            enqueueSafeMergeClausePropagees(reason[q],
                                            updated == v ? w :v,
                                            getAssignment(updated));
          }
          continue;
        }

        Lit falselit = (res == l_True) ? Lit(q, true) : Lit(q, false) ;
        vec<Clause*>& ws = watches[toInt(~falselit)];
        for (i = j = static_cast<Clause**>(ws), end = i + ws.size();  i != end;) {
          Clause& c = **i++;
          // Make sure the false literal is data[1]:
          if (c[0] == falselit) {
            c[0] = c[1];
            c[1] = falselit;
          }
          assert(c[1] == falselit);
          // If 0th watch is true, then clause is already satisfied.
          Lit first = c[0];
          factorMerge(first, v, w);
          if (valueOfMerge(v, w, sign(first)) == l_True) {
            *j++ = &c;
          } else {
            // Look for new watch:
            for (int k = 2; k < c.size(); k++) {
              if (canBeUsedAsWatch(c[k])) {
                c[1] = c[k];
                c[k] = falselit;
                watches[toInt(~c[1])].push(&c);
                Var q_ = var(c[1]);
                factorMerge(q_, v_, w_);
                if (locations[q_].first == NULL) {
                  primaryPointers[v_].push(q_);
                  locations[q_].first = &(primaryPointers[v_].last());
                }
                if (locations[q_].second == NULL) {
                  primaryPointers[w_].push(q_);
                  locations[q_].second = &(primaryPointers[w_].last());
                }
                goto FoundVertexWatch;
              }
            }

            // Did not find watch -- clause is unit under assignment:
            *j++ = &c;
            if (valueOfMerge(first) == l_False) {
              confl = &c;
              qhead = trail.size();
              // Copy the remaining watches:
              while (i < end) {
                *j++ = *i++;
              }
             } else {
              if (isAssigned(v)) {
                enqueueSafeMergeClausePropagees(&c, w, getAssignment(v));
              } else if (isAssigned(w)) {
                enqueueSafeMergeClausePropagees(&c, v, getAssignment(w));
              } else if (value(c[0]) == l_Undef){
                uncheckedEnqueue(c[0], &c);
              }
            }
          }
          FoundVertexWatch: {}
        }
        ws.shrink(i - j);
        factorMerge(q, v, w);
        if ((watches[toInt(Lit(q, true))].size() + watches[toInt(Lit(q, false))].size()) == 0) {
          bool updatedIsV = (v == updated);
          if ( updatedIsV ) {
            *locations[q].second = primaryPointers[w].last();
            Var q_ = *locations[q].second;
            factorMerge(q_, v_, w_);
            if (v_ == w) {
              locations[q_].first = locations[q].second;
            } else {
              locations[q_].second = locations[q].second;
            }
            primaryPointers[w].pop();
          } else {
            *locations[q].first =  primaryPointers[v].last();
            Var q_ = *locations[q].first;
            factorMerge(q_, v_, w_);
            if (v_ == v) {
              locations[q_].first = locations[q].first;
            } else {
              locations[q_].second = locations[q].first;
            }
            primaryPointers[v].pop();
          }
          locations[q] = std::make_pair<Var*, Var*>(NULL, NULL);
        } else if (locations[q].first == NULL) {
          primaryPointers[v].push(q);
          locations[q].first = &(primaryPointers[v].last());
        } else if (locations[q].second == NULL) {
          primaryPointers[w].push(q);
          locations[q].second = &(primaryPointers[w].last());
        }
      }
    }

    for (int x = 0; x < vertexStack.size(); x++) {
      vertexSeen[vertexStack[x]] = false;
    }
    vertexStack.clear();
  } while (qhead < trail.size());

  propagations += num_props;
  simpDB_props -= num_props;
  return confl;
}
/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt {
  bool operator () (Clause* x, Clause* y) {
    return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity());
  }
};
void Solver::reduceDB()
{
    assert(decisionLevel() == 0);
    int     i, j;
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

    sort(learnts, reduceDB_lt());
    for (i = j = 0; i < learnts.size() / 2; i++){
        if (learnts[i]->size() > 2 && !locked(*learnts[i])) {
            removeClause(*learnts[i]);
        } else {
            learnts[j++] = learnts[i];
        }
    }
    for (; i < learnts.size(); i++){
        if (learnts[i]->size() > 2 && !locked(*learnts[i]) && learnts[i]->activity() < extra_lim) {
            removeClause(*learnts[i]);
        } else {
            learnts[j++] = learnts[i];
        }
    }
    learnts.shrink(i - j);
}

void Solver::removeSatisfied(vec<Clause*>& cs)
{
    int i,j;
    for (i = j = 0; i < cs.size(); i++){
        if (satisfied(*cs[i]))
            removeClause(*cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != NULL)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    //removeSatisfied(learnts);
//    if (remove_satisfied)        // Can be turned off.
//        removeSatisfied(clauses);

    // Remove fixed variables from the variable heap:
    order_heap.filter(VarFilter(*this));

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]
|
|  Description:
|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
|    indicate infinity.
|
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts, int nof_learnts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;

    starts++;
    bool first = true;

    unsigned int colorInfo = 0;
    bool seenSigned = false;
    Vertex freeVertex = 0;
    for (;;){
        Clause* confl = propagate();
        if (confl != NULL){
            // CONFLICT
            conflicts++; conflictC++; first = false;
            if (decisionLevel() == 0) {
              return l_False;
            }

            bool first = true;
            int res = vertex_undef;
            do {
              learnt_clause.clear();
              analyze(confl, learnt_clause, backtrack_level);
              res = transformConflict(learnt_clause,
                                      backtrack_level,
                                      freeVertex,
                                      seenSigned,
                                      colorInfo);
              if (res == vertex_unsat) {
                return l_False;
              } else if (res == vertex_reitarate) {
                cancelUntil(backtrack_level);
                if (!first) {
                  free(confl);
                }
              }
              if (learnt_clause.size() > 1)
                  confl = Clause_new(learnt_clause, true);
              first = false;
            } while (res == vertex_reitarate);

            if (learnt_clause.size() == 0)
                return l_Undef;
            else if (learnt_clause.size() == 1) {
                cancelUntil(0);
                if (value(learnt_clause[0]) == l_Undef)
                    uncheckedEnqueue(learnt_clause[0]);
            } else {
                Clause* c = confl;

                attachMergeClause(c);
                cancelUntil(backtrack_level);
                learnts.push(c);
#ifdef SHOWLEARNT
                reportf("\nLearnt:%d\n", (int)c);
                printClauseFancy(*c);
                reportf("\n\n");
#endif
                enqueueUnsafeMergeClausePropagees(
                    c, freeVertex, seenSigned, colorInfo);
            }
            vertexDecayActivity();
            claDecayActivity();

        }else{
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (nof_learnts >= 0 && learnts.size()-nAssigns() >= nof_learnts) {
                // Reduce the set of learnt clauses:
                cancelUntil(0);
                reduceDB();
            }

            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                next = pickBranchLit(polarity_mode, random_var_freq);

                if (next == lit_Undef)
                    // Model found:
                    return l_True;
            }

            // Increase decision level and enqueue 'next'
            assert(value(next) == l_Undef);
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}


double Solver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}


bool Solver::solve(const vec<Lit>& assumps)
{
    model.clear();
    conflict.clear();

    if (!ok) return false;

    assumps.copyTo(assumptions);

    double  nof_conflicts = restart_first;
    double  nof_learnts   = nClauses() * learntsize_factor;
    lbool   status        = l_Undef;

    if (verbosity >= 1){
        reportf("============================[ Search Statistics ]==============================\n");
        reportf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        reportf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        reportf("===============================================================================\n");
    }

    // Search:
    while (status == l_Undef){
        if (verbosity >= 1)
            reportf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", (int)conflicts, order_heap.size(), nClauses(), (int)clauses_literals, (int)nof_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progress_estimate*100), fflush(stdout);
        status = search((int)nof_conflicts, (int)nof_learnts);
        nof_conflicts *= restart_inc;
        nof_learnts   *= learntsize_inc;
    }

    if (verbosity >= 1)
        reportf("===============================================================================\n");


    if (status == l_True){
        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
#ifndef NDEBUG
        verifyModel();
#endif
    }else{
        assert(status == l_False);
        if (conflict.size() == 0)
            ok = false;
    }

    cancelUntil(0);
    return status == l_True;
}

//=================================================================================================
// Debug methods:


void Solver::verifyModel()
{
    bool failed = false;
    for (int i = 0; i < clauses.size(); i++){
        assert(clauses[i]->mark() == 0);
        Clause& c = *clauses[i];
        for (int j = 0; j < c.size(); j++)
            if (modelValue(c[j]) == l_True)
                goto next;

        reportf("unsatisfied clause: ");
        printClause(*clauses[i]);
        reportf("\n");
        failed = true;
    next:;
    }

    assert(!failed);

    // reportf("Verified %d original clauses.\n", clauses.size());
}


void Solver::checkLiteralCount()
{
    // Check that sizes are calculated correctly:
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (clauses[i]->mark() == 0)
            cnt += clauses[i]->size();

    if ((int)clauses_literals != cnt){
        fprintf(stderr, "literal count: %d, real value = %d\n", (int)clauses_literals, cnt);
        assert((int)clauses_literals == cnt);
    }
}

void Solver::dumpColorAssignments() {
  reportf("dumping color assingments, at level:%d ..\n", decisionLevel());
  for (int i =0; i < numberOfVertices; i++) {
    if (isAssigned(i)) {
      continue;
    }
    reportf("vertex:%03d -> c:%03d ..\n", i, colorAssignedMask[i]);
  }
}


void Solver::dumpTrail(Clause * confl, vec<bool>& mySeen) {  // NOLINT
  reportf("\n");
  reportf("dumping conflict:%lld at level:%d\n", conflicts, decisionLevel());
  for (int i =0; i < trail.size(); i++) {
    if (((!isVertex(trail[i]) && (confl != NULL)) || !mySeen[var(trail[i])]) && level[var(trail[i])] != decisionLevel()) {
      continue;
    }
    reportf("[%d] ", i);
    printLitFancy(trail[i]);
    reportf("\t\ts[%d]\t >",
            seen[var(trail[i])]);
    if (reason[var(trail[i])] != NULL) {
      Clause& clause = *reason[var(trail[i])];
      printClauseFancy(clause);
    } else {
      reportf("No reason...");
    }
    reportf("\n");
  }
  reportf("\n\n");

  if (confl != NULL) {
    reportf("\n\nResulted in the falsefied clause:%d at:%d ..\n", confl->learnt(), decisionLevel());
    printClauseFancy(*confl);
    reportf("\n\n");
  }

}
