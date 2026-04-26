#include "internal.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>


namespace CaDiCaL {

  /*
    Add a clause to the formula but mark it as garbage.
  
    This prevents it from interacting with the following encoding steps (we skip garbage clauses).

    Some of these clauses will be kept (marked as not garbage), 
    if reencodeKeep is set to true e.g., (yi -> yi+1),
    or elim is set to false.
  */
Clause * Internal::add_proof_clause (const vector<int> &lits, Clause *temp, bool garbage) {
  for (const auto &lit : lits)
    clause.push_back (lit);
  
  Clause *c = new_clause_as (temp);
  clause.clear ();
  c->garbage = garbage;
  return c;
}

/*
  Add the linear encoding for a ULC so that we can 
  add the sequential counter. These clauses will be 
  deleted after the sequential counter is done. 
  They need to be deleted in reverse order to be RAT 
  steps, and will be added to the extension stack.
*/
void Internal::add_linear_encoding (vector<int> lits, vector<Clause *> &clauses, Clause * c) {

  // base case, less than 5 literals
  if (lits.size () < 5) {
    // pairwise encoding 
    for (unsigned i = 0; i < lits.size () - 1; i++) {
      for (unsigned j = i + 1; j < lits.size (); j++) {
        clauses.push_back (add_proof_clause (vector<int>{-lits[j], -lits[i]}, c, true));
      }
    }
    return;
  }

  // recursive case
  int aux = max_var + 1;
  external->internalize (aux);

  // paiwise encoding for first three literals 
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = i + 1; j < 3; j++) {
      clauses.push_back (add_proof_clause (vector<int>{-lits[j], -lits[i]}, c, true));
    }
    clauses.push_back (add_proof_clause (vector<int>{-aux, -lits[i]}, c, true));
  }

  // clause for backwards propagation of aux
  clauses.push_back (add_proof_clause (vector<int>{aux, lits[0], lits[1], lits[2]}, c, true));

  // recursive call with first three literals removed and aux added to front of lits
  vector<int> new_lits;
  new_lits.push_back (-aux);
  new_lits.insert (new_lits.end (), lits.begin () + 3, lits.end ());
  add_linear_encoding (new_lits, clauses, c);

}

/*
  Delete the linear encoding from above, taking care with RAT steps
  Use old_max_var to determine which clauses contain new variables
  added to the linear encoding.
*/
void Internal::delete_linear_encoding (vector<Clause *> &clauses, int old_max_var) {

  vector<int> elim_vars;

  // delete clauses in reverse order
  for (int i = clauses.size () - 1; i >= 0; i--) {
    Clause * c = clauses[i];
    c->garbage = false;

    // sort clause largest to smallest by absolute value
    sort (c->literals, c->literals + c->size, [&](int a, int b) {
      return abs(a) > abs(b);
    });

    if (abs(c->literals[0]) > old_max_var) {
      // special case for RAT clauses
      external->push_clause_on_extension_stack (c, c->literals[0]);
      elim_vars.push_back (abs (c->literals[0]));
    }

    mark_garbage (c);
  }

  for (auto v : elim_vars) {
    if (flags (v).eliminated ()) continue;
    mark_eliminated (v);
  }
}

/*
  Add remaining edges of the pairiwse encoding
  Do not add duplicate clauses, we need to check this 
  because of hybrid XLCs. Add to coresponding maps 
  for quick lookup and deletion.
*/
void Internal::add_pairwise_encoding (vector<int> lits, vector<set<int>> &bin_edges, vector<set<Clause*>> &bin_cls, Clause * c) {
  if (lits.size() <= 1) return;
  for (unsigned i = 0; i < lits.size () - 1; i++) {
    for (unsigned j = i + 1; j < lits.size (); j++) {
      int l1 = lits[i];
      int l2 = lits[j];

      // check if clause already exists
      if (bin_edges[vlit(-l1)].find(-l2) != bin_edges[vlit(-l1)].end()) continue;

      // for RAT check on unqiue literal
      if (noccs((l2)) == 1) swap (l1, l2);

      Clause * cN = add_proof_clause (vector<int>{-l1, -l2}, c, false);
      bin_cls[vlit(-l1)].insert (cN);
      bin_cls[vlit(-l2)].insert (cN);
      bin_edges[vlit(-l1)].insert(-l2); 
      bin_edges[vlit(-l2)].insert(-l1);
    }
  }
}

/*
  Delete the pairwise encoding from above, all are RUP 
  because of the sequential counter.
*/
void Internal::delete_pairwise_encoding (vector<int> lits, vector<set<Clause*>> &bin_edges_vec_cls) {
  for (auto lit : lits) {
    for (auto lit2 : lits) {
      if (lit == lit2) continue;
      // find binary clause -lit, -lit2
      for (auto c : bin_edges_vec_cls[vlit(-lit)]) {
        if (c->literals[0] == -lit2 || c->literals[1] == -lit2) {
          if (!c->garbage)
            mark_garbage (c);
          break;
        }
      }
    }
  }
}

/*
  Read in the cliques from a file and store them in a vector of vectors
  where each vector is a clique.

  File name is hard coded to cliques.txt. Each line in the file is a clique with the literals separated by spaces.

  Unused currently
*/
void read_cliques (unordered_map<int,int> &cliques) {
  ifstream file ("cliques.txt");
  string line;
  int clique_num = 0;
  while (getline (file, line)) {
    istringstream iss (line);
    int n;
    while (iss >> n) {
      // if n is in cliques skip
      if (cliques.find(n) == cliques.end()) {
        cliques[n] = clique_num;
      }
    }
    clique_num++; 
  }
  // print out parsed cliques
  // for (auto const& x : cliques)
  // {
  //   cout << x.first  // string (key)
  //        << ':'
  //        << x.second // string's value 
  //        << std::endl ;
  // }
}

/*
  Alignment of xlcs using binary clauses.
*/
void Internal::align_ulc_clauses (unordered_map<uint64_t,int> &ulcs, vector<set<int>>  &bin_edges) {

  vector<int> var_alignment;
  var_alignment.resize (2 * (1 + vsize), -1);
  
  // Get xlc clauses into a vector
  vector<Clause*> ulc_clauses;
  for (auto c : clauses) {
    if (c->garbage) continue;
    if (ulcs.find(c->id) != ulcs.end()) {
      for (auto lit: *c) { // remove from bin map edges from same ULC
        // do not align ulcs that fail to meet the cutoff
        // since we will not be reencoding them
        
        var_alignment [vlit(lit)] = 0; // mark variables needing alignment
        for (auto other_lit : *c) {
          if (lit == other_lit) continue;
          if (bin_edges[vlit(-lit)].find(-other_lit) != bin_edges[vlit(-lit)].end()) {
            bin_edges[vlit(-lit)].erase(-other_lit);
            bin_edges[vlit(-other_lit)].erase(-lit);
          }
        }
      }
      if (c->size <= opts.orderencodecut) continue;
      ulc_clauses.push_back(c);
    }
  }

  // sort ulc clauses by size, largest to smallest
  sort (ulc_clauses.begin(), ulc_clauses.end(), [](Clause* a, Clause* b) {
    return a->size > b->size;
  });

  cout << "c Aligning " << ulc_clauses.size () << " clauses" << endl;

  int ulc_alignment = 0; // alignment type (count of literals aligned)
  int alignment = 0;     // alignment column

  // loop over all ulc clauses and align literals
  for (auto c : ulc_clauses) {
    // go over each literal in the clause, checking if it is already aligned, otherwise add it to the next alignment

    // sort literals in c (so natural is default outcome)
    sort (c->literals, c->literals + c->size, [](int a, int b) {
      return abs(a) < abs(b); // sort by absolute value
    });

    for (auto lit : *c) {
      if (var_alignment[vlit (lit)] > 0) continue; // already aligned
      alignment++;
      var_alignment[vlit (lit)] = alignment;

      vector<int> newly_aligned;
      newly_aligned.push_back (lit);
      while (newly_aligned.size()) { // go through bin map and align all connected lits
        int new_lit = newly_aligned.back();
        newly_aligned.pop_back();
        for (auto other_lit : bin_edges[vlit (-new_lit)]) {
          if (var_alignment[vlit (-other_lit)] == 0) {
            var_alignment[vlit (-other_lit)] = alignment;
            newly_aligned.push_back (-other_lit);
            if (ulc_alignment >= 0)
              ulc_alignment++;
          } 
        }
      }
    }

    // Everything now aligned and we can sort the lits in the clause
    // lexicographic ordering with abs(var) if alignment is the same
    sort (c->literals, c->literals + c->size, [&](int a, int b) {
      // CHANGE!
      if (var_alignment[vlit (a)] == var_alignment[vlit (b)])
        return abs(a) < abs(b); // if both are unaligned sort by abs value
      return var_alignment[vlit (a)] < var_alignment[vlit (b)];
    });

    // check if any literals share same alignment (unaligneable in this case)
    for (int i = 0; i < c->size - 1; i++) {
      int lit = c->literals[i];
      int next_lit = c->literals[i+1];
      if (var_alignment[vlit (lit)] == var_alignment[vlit (next_lit)])
        ulc_alignment = -1;
    }

  }

  stats.order.ulc_alignment = ulc_alignment; // -1 unalignable, 0 no alignment, >0 aligned
}

void Internal::mark_unique_literal_clauses (unordered_map<uint64_t,int> &ulcs, vector<set<int>>  &bin_edges) {

  // with the options we will decide how to categorize
  int ordertype = opts.orderencodetype;
  // 1: ulc, 2: all bin_edges, 3: xlc:ulc,bin_edges,hybrid

  for (auto c : clauses) {
    if (c->garbage) continue;

    int ul_cnt = 0;
    int bin_edge_cnt = 0;
    bool can_reencode = false;
    bool is_pure = false;

    vector<int> bin_lits;

    if (ordertype != 2) {
      can_reencode = true; // default to true, we will check if we can reencode
      for (auto lit : *c) {
        if (!noccs (-lit)) { // no reason to reencode pure literal (satisfies clause...)
          is_pure = true;
          can_reencode = false;
          break;
        }
        if (noccs (lit) != 1) {
          can_reencode = false;
          if (ordertype == 1) {
            break;
          } else { // could by hybrid
            bin_lits.push_back (lit);
          }
        } else ul_cnt++;
      }
    }

    if (is_pure) continue;

    // bin_edge check necessary in all cases
    bool is_bin_edge = true;
    for (int i = 0; i < c->size - 1; i++) {
      for (int j = i + 1; j < c->size; j++) {
        if (bin_edges[vlit(-c->literals[i])].find(-c->literals[j]) == bin_edges[vlit(-c->literals[i])].end()) {
          is_bin_edge = false;
          break;
        }
      }
    }
    if (is_bin_edge) {
      bin_edge_cnt = c->size;
      can_reencode = true;
    } else if (!can_reencode && ordertype == 3) {
      // perform hybrid check 
      is_bin_edge = true;
      for (unsigned i = 0; i < bin_lits.size () - 1; i++) {
        for (unsigned j = i + 1; j < bin_lits.size (); j++) {
          if (bin_edges[vlit(-bin_lits[i])].find(-bin_lits[j]) == bin_edges[vlit(-bin_lits[i])].end()) {
            is_bin_edge = false;
            break;
          }
        }
      }
      if (is_bin_edge) {
        bin_edge_cnt = bin_lits.size ();
        can_reencode = true;
      }
    }

    bool is_ulc = (ul_cnt == c->size); // all literals are ulc
    bool is_bin = (bin_edge_cnt == c->size); // all literals are bin_edges

    if (can_reencode) {
      int re_type = -1;
      if (is_ulc && bin_edge_cnt == 0) {
        re_type = 1;
        if (c->size > opts.orderencodecut)
          stats.order.pre_ulcs++;
        else
          stats.order.pre_ulcs_below++;
      } else if (is_ulc && is_bin) {
        re_type = 4;
        if (c->size > opts.orderencodecut)
          stats.order.pre_both++;
        else 
          stats.order.pre_both_below++;
      } else if (!ul_cnt) {
        re_type = 2;
        if (c->size > opts.orderencodecut)
          stats.order.pre_bin++;
        else 
          stats.order.pre_bin_below++;
      } else if (ul_cnt && is_bin) {
        re_type = 5;
        if (c->size > opts.orderencodecut)
          stats.order.pre_hybf++;
        else 
          stats.order.pre_hybf_below++;
      } else {
        re_type = 3;  
        if (c->size > opts.orderencodecut)
          stats.order.pre_hyb++;
        else 
          stats.order.pre_hyb_below++;
      }

      if (ordertype == 1) {
        if (re_type != 1 && re_type != 4)
          continue;
      }
      ulcs[c->id] = re_type;
    }
  }
}


// Not Used
void Internal:: xlc_resolutions (unordered_map<uint64_t,int> &ulcs) {
  
  vector<vector<Clause*>> xlc_occ, var_occ;
  vector<int> stamp, in_xlc;
  stamp.resize ((vsize + 1) * 2, 0);
  in_xlc.resize ((vsize + 1) * 2, 0);
  int time = 1;

  LOG ("Get clause pointers");
  xlc_occ.resize ((vsize + 1) * 2, vector<Clause*>());
  var_occ.resize ((vsize + 1) * 2, vector<Clause*>());
  
  // loop over all clauses
  for (auto c : clauses) {
    if (c->garbage) continue;
    // if xlc
    if (ulcs.find(c->id) != ulcs.end()) {
      for (auto lit : *c) {
        xlc_occ [vlit (lit)].push_back (c);
        in_xlc [vlit (lit)] = 1;
      }
    } else {
      for (auto lit : *c) {
        var_occ [vlit (lit)].push_back (c);
      }
    }
  }

  vector<bool> forced_vars;
  forced_vars.resize ((vsize + 1), 0);
  for (auto idx: vars) {
    if (xlc_occ[vlit(idx)].size() == 1 && xlc_occ[vlit(-idx)].size() == 1 && var_occ[vlit(idx)].size() == 0 && var_occ[vlit(-idx)].size() == 0) {
      forced_vars[abs(idx)] = true; // this variable is involved in an XLC
    }
  }

  // Get variables that are forced into resolution
  LOG ("Get forced resolutions");
  for (auto idx: vars) {
    if (val (idx)) continue;

    if (in_xlc [vlit (idx)] == 0 || in_xlc [vlit (-idx)] == 0) continue; // does not occur in any xlc
    
    bool pos_c = false, neg_c = false;

    if (forced_vars [idx]) {
      pos_c = neg_c = true;
    }

    // cout << "checking Variable " << idx << " for forced resolution" << endl;
    // check if clash, if so perform resolutions (not elimination if not purely ulc)

    //
    // for (auto c : xlc_occ [vlit (idx)]) {
    //   if (c->garbage) continue;
    //   pos_c = true;
    //   break;
    // }
    // for (auto c : xlc_occ [vlit (-idx)]) {
    //   if (c->garbage) continue;
    //   neg_c = true;
    //   break;
    // }

    

    // check if any positive clauses ulc, and any negative clauses ulc
    if (pos_c && neg_c) {

      // cout << "c Eliminating clashing XLCs for variable " << idx << endl;

      stats.order.clashElim++;

      vector<Clause*> pos_clss, neg_clss;
      // bool pos_pres = false, neg_pres = false;
      // if (opts.orderencodetype == 2 || opts.orderencodetype == 3) {
      //   for (auto c : var_occ [vlit (idx)]) {
      //     if (c->garbage) continue;
      //     pos_pres = true;
          
      //       pos_clss.push_back (c);
      //       break;
      //   }
      //   for (auto c : var_occ [vlit (-idx)]) {
      //     if (c->garbage) continue;
      //     neg_pres = true;
          
      //       neg_clss.push_back (c);
      //       break;
      //   }
      // } 
      // if (pos_clss.size() > 0 || neg_clss.size () > 0) 
      //   continue; // early exit

      // cout << "c var " << idx << endl;
      // cout << "c pos " << pos_clss.size () << " neg " << neg_clss.size () << endl;
      // cout << "c pos XLCs " << xlc_occ [vlit (idx)].size () << " neg XLCs " << xlc_occ [vlit (-idx)].size () << endl;

      // perform resolution on clashing XLCs
      // resolvent is an XLC
      for (auto pos_c : xlc_occ [vlit (idx)]) {
        if (pos_c->garbage) continue;
        int ptype = ulcs[pos_c->id];
        for (auto neg_c : xlc_occ [vlit (-idx)]) {
          if (pos_c->garbage || neg_c->garbage) continue;
          int ntype = ulcs[neg_c->id];
          stats.order.nres++;

          // cout << "c Resolving XLC " << pos_c->id << " (type " << ptype << ") with XLC " 
          //      << neg_c->id << " (type " << ntype << ")" << endl;

          // resolve pos_c and neg_c; 
          bool is_taut = false;
          int restype = 0;
          time += 1;

          for (auto lit : *pos_c) {
            if (abs(lit) == idx) continue;
            if (flags (abs(lit)).eliminated ())
              continue;
            const signed char tmp = val (lit);
            if (tmp > 0) {
              is_taut = lit;
              break;
            } else if (tmp < 0) continue;
            clause.push_back (lit);
            stamp [vlit (lit)] = time;
          }
          for (auto lit : *neg_c) {
            if (abs(lit) == idx) continue;
            if (flags (abs(lit)).eliminated ())
              continue;
            const signed char tmp = val (lit);
            if (tmp > 0) {
              is_taut = lit;
              break;
            } else if (tmp < 0) continue;
            if (stamp [vlit (-lit)] == time) {
              is_taut = true;
              break;
            }
            else if (stamp [vlit (lit)] == time)
              continue;
            clause.push_back (lit);
            stamp[vlit (lit)] = time;
          }

          if (is_taut) {
            stats.order.tautres++;  
          } else if (clause.size() > 1) {
            Clause * res = new_clause_as (pos_c);
            // Add to XLCs new clauses
      // Add to map new XLCs
            if (ptype == ntype)
              restype = ptype;
            else
              restype = 3;
            if (restype == 4) restype = 1;
            ulcs[res->id] = restype;
            for (auto lit : clause) {

              // only works for ulc clash
              if (xlc_occ [vlit (lit)].size () == 1) {
                xlc_occ [vlit (lit)][0] = res;
              } else 
                xlc_occ [vlit (lit)].push_back (res);
            }
          } else if (clause.size() == 1) {
            assign_unit (clause[0]);
            if (!propagate ()) {
              learn_empty_clause ();
              return;
            }
        //   LOG ("Unit clause %d found", clause[0]);
            stats.order.unitsres++;
          } else {
            conflict = pos_c;
            learn_empty_clause ();
            conflict = 0;
            return;
          }
          clause.clear ();
        }
      }

      // need to perform resolutions on these clashing clauses
      for (auto pos_c : pos_clss) {
        if (pos_c->garbage) continue;
        for (auto neg_c : neg_clss) {
          if (pos_c->garbage || neg_c->garbage) continue;
          stats.order.nNonXlcRes++;
          // resolve pos_c and neg_c; 
          bool is_taut = false;
          time += 1;

          for (auto lit : *pos_c) {
            if (abs(lit) == idx) continue;
            if (flags (abs(lit)).eliminated ())
              continue;
            const signed char tmp = val (lit);
            if (tmp > 0) {
              is_taut = lit;
              break;
            } else if (tmp < 0) continue;
            clause.push_back (lit);
            stamp [vlit (lit)] = time;
          }
          for (auto lit : *neg_c) {
            if (abs(lit) == idx) continue;
            if (flags (abs(lit)).eliminated ())
              continue;
            const signed char tmp = val (lit);
            if (tmp > 0) {
              is_taut = lit;
              break;
            } else if (tmp < 0) continue;
            if (stamp [vlit (-lit)] == time) {
              is_taut = true;
              break;
            }
            else if (stamp [vlit (lit)] == time)
              continue;
            clause.push_back (lit);
            stamp[vlit (lit)] = time;
          }

          if (is_taut) {
            stats.order.tautres++;  
          } else if (clause.size() > 1) {
            Clause * res = new_clause_as (pos_c);
            for (auto lit : clause) {
              var_occ [vlit (lit)].push_back (res);
            }
          } else if (clause.size() == 1) {
            assign_unit (clause[0]);
            if (!propagate ()) {
              learn_empty_clause ();
              return;
            }
        //   LOG ("Unit clause %d found", clause[0]);
            stats.order.unitsres++;
          } else {
            conflict = pos_c;
            learn_empty_clause ();
            conflict = 0;
            return;
          }
          clause.clear ();
        }
      }

      // // cout << "c pre erasures" << endl;
      // Turn the vector occ into aisngle occ
      // eliminating way more variables?
      
      // Erase from XLCs old clauses
      for (auto c : xlc_occ [vlit (idx)]) {
        if (c->garbage) continue;
        pos_clss.push_back (c);
      }
      for (auto c : xlc_occ [vlit (-idx)]) {
        if (c->garbage) continue;
        neg_clss.push_back (c);
      }

      // delete with pivot for extension stack
      for (auto pos_c : pos_clss) {
        if (pos_c->garbage) continue;
        external->push_clause_on_extension_stack (pos_c, idx);
        mark_garbage (pos_c);
      }
      for (auto neg_c : neg_clss) {
        if (neg_c->garbage) continue;
        external->push_clause_on_extension_stack (neg_c, -idx);
        mark_garbage (neg_c);
      }
      mark_eliminated (abs(idx));

      // cout << "c Eliminated clashing XLCs for variable " << idx << endl;

    }
  }

}

// Resolutions (two passes over clauses)
  //  1st, mark all variables with single pos and neg occurence
  //  2nd, resolve all ulcs clauses with these variables
void Internal:: order_encode_resolutions (unordered_map<uint64_t,int> &ulcs) {
  
  vector<int> proc_stack, forced, uniq;
  vector<int> res_pivot_list;
  vector<Clause*> res_list;
  vector<Clause*> var_occ;
  
  vector<int> stamp;
  stamp.resize ((vsize + 1) * 2, 0);
  int time = 1;

  LOG ("Get clause pointers");
  var_occ.resize ((vsize + 1) * 2, 0);
  for (auto c : clauses) {
    if (c->garbage) continue;
    for (auto lit : *c) {
      var_occ [vlit (lit)] = c;
    }
  }

  // Get variables that are forced into resolution
  LOG ("Get forced resolutions");
  forced.resize ( (vsize + 1), 0);
  uniq.resize ( (vsize + 1) * 2, 0);
  for (auto idx: vars) {
    if (noccs (idx) == 1 && noccs (-idx) == 1) {
      forced[idx] = 1;
      LOG ("Found forced resolution for variable %d",idx);
      stats.order.forced_res++;
    }
  }

  auto old_id = clauses.back()->id;

  for (unsigned i = 0; i < clauses.size(); i++) {
    Clause *c = clauses[i];
    if (c->garbage) continue;
    if (ulcs.find(c->id) == ulcs.end()) continue;
    if (ulcs[c->id] != 1 && ulcs[c->id] != 4) continue; // only ulc and Both

    int forced_var = 0;
    for (auto lit: *c) {
      if (forced [abs(lit)] == 1 ) {
        if (flags (abs(lit)).eliminated ())
              continue;
        if (var_occ [vlit (-lit)]->garbage) {
          // now a pure literal
          external->push_clause_on_extension_stack (c, lit);
          mark_garbage (c);
          mark_eliminated (abs(lit));
          if (c->id <= old_id) // counted here as technically a resolution even if result of tautology earlier
            stats.order.res_old++;
          else 
            stats.order.res_new++;
          stats.order.nres++;
          break;
        }
        if (ulcs.find(var_occ [vlit (-lit)]->id) != ulcs.end() && (ulcs[var_occ [vlit (-lit)]->id] == 1 || ulcs[var_occ [vlit (-lit)]->id] == 4)) {
          forced_var = abs(lit);
          break;
        }
      }
    }

    if (forced_var) {
      // resolution candidate found
      bool is_taut = false;
      int restype = 0;
      time += 1;
      forced_var = 0;
      for (auto lit: *c) {
        if (flags (abs(lit)).eliminated ())
              continue;
        if (!forced_var && forced [abs(lit)] == 1 ) {
          if (!var_occ [vlit (-lit)]->garbage && ulcs.find(var_occ [vlit (-lit)]->id) != ulcs.end() && (ulcs[var_occ [vlit (-lit)]->id] == 1 || ulcs[var_occ [vlit (-lit)]->id] == 4)) {
            forced_var = lit;
            proc_stack.push_back (-lit);
          }
        }
        clause.push_back (lit);
        stamp [vlit (lit)] = time;
      }

      restype = ulcs[c->id];
      res_list.push_back (c);
      res_pivot_list.push_back (forced_var);
      // c->garbage = true;

      if (c->id <= old_id)
        stats.order.res_old++;
      else 
        stats.order.res_new++;

      while (proc_stack.size () > 0) {

        int next_lit = proc_stack.back ();
        proc_stack.pop_back ();

        LOG ("Popping %d from processed stack",next_lit);

        Clause *next_c = var_occ [vlit (next_lit)];

        if (next_c->id <= old_id)
          stats.order.res_old++;
        else 
          stats.order.res_new++;

        // remove -next_lit from clause
        clause.erase (find (clause.begin(), clause.end(), -next_lit));

        // perform resolution
        for (auto lit : *next_c) {
          if (flags (abs(lit)).eliminated ())
              continue;
          if (abs(lit) == abs(next_lit)) continue;
          if (stamp [vlit (lit)] == time)
            continue;
          else if (stamp [vlit (-lit)] == time)
            is_taut = true;
          // if (forced [abs(lit)] == 1) {
          //   if (ulcs.find(var_occ [vlit (-lit)]->id) != ulcs.end() && !var_occ [vlit (-lit)]->garbage) {
          //     proc_stack.push_back (-lit);
          //     LOG ("Pushing %d to processed stack", abs(lit));
          //     // cout << "Pushing " << (-lit) << " to processed stack" << endl;
          //   }
          // }
          clause.push_back (lit);
          stamp[vlit (lit)] = time;
        }

        restype = 1;
        stats.order.nres++;
        // next_c->garbage = true;
        res_list.push_back (next_c);
        res_pivot_list.push_back (next_lit);

        break;
      }

      // cout << "c Resolvent: ";
      // for (auto lit: clause) {
      //   cout << lit << " ";
      
      // } cout << endl;

      assert (clause.size () || is_taut);
      LOG (clause, "New resolvent");
      if (!is_taut) {
        assert (clause.size () > 1);
        // c->garbage = false;
        Clause * res = new_clause_as (c);
        ulcs[res->id] = restype;
        for (auto lit : clause) {
          var_occ [vlit (lit)] = res;
        }
      } else 
        stats.order.tautres++;

      // for (unsigned i = 0; i < res_list.size(); i++)
      //   res_list[i]->garbage = false;

      // delete all clauses that were resolved
      for (unsigned i = 0; i < res_list.size(); i++) {
        int pivot = res_pivot_list[i];
        Clause *res_clause = res_list[i];
          if (res_clause->garbage) continue;
        // if (proof)
        //   proof->weaken_minus (res_clause); // why weaken minus?
        external->push_clause_on_extension_stack (res_clause, pivot);
        mark_garbage (res_clause);
      }

      for (unsigned i = 0; i < res_list.size(); i++) {
        int pivot = res_pivot_list[i];
        if (!flags (abs(pivot)).eliminated ())
          mark_eliminated (abs(pivot));
      }
      
      proc_stack.clear();
      res_list.clear();
      res_pivot_list.clear();
      clause.clear ();
    }
  }
}

/*
  Perform order encoding on the formula.

  The order encoding is performed in two steps:
    1. Resolutions
    2. Order Encoding

  This is broken down further as:
    1. Resolutions
      1.1. Get occurence counts
      1.2. Get forced resolutions
      1.3. Perform resolutions on forced variables
    2. Order Encoding
      2.1. Get occurence counts
      2.2. Find substitutions for order encoding
      2.3. Remove tautologies
      2.4. Reencode clauses
    3. Delete proof only clauses
    4. Delete garbage clauses
    5. Reset occurence counts
    6. Clear and reconnect watches
*/
// hash function for pairs


void Internal::order_encode () {
  START (order);
  LOG ("START order encoding");

  // mark_satisfied_clauses_as_garbage(); // in case we need to remove some units
  // reduce ();
  garbage_collection ();


  vector<set<int>> bin_edges_vec;
  vector<set<Clause*>> bin_edges_vec_cls;
  
  bin_edges_vec.resize (2 * (vsize + 1));
  bin_edges_vec_cls.resize (2 * (vsize + 1));
  for (auto c : clauses) {
    if (c->garbage) continue;
    if (c->size == 2) {
      bin_edges_vec[vlit(c->literals[0])].insert((c->literals[1]));
      bin_edges_vec[vlit(c->literals[1])].insert((c->literals[0]));
      bin_edges_vec_cls[vlit(c->literals[0])].insert(c);
      bin_edges_vec_cls[vlit(c->literals[1])].insert(c);
    }
  }

  cout << "c Made edge map" << endl;

  LOG ("Get occurence counts");
  init_noccs ();
  for (auto c : clauses) {
    if (c->garbage) continue;
    for (auto lit : *c) {
      noccs (lit)++;
    }
  }

  vector<int> stamp;
  stamp.resize ((vsize + 1) * 2, 0);
  int time = 1;

  unordered_map<uint64_t,int> ulcs;
  mark_unique_literal_clauses (ulcs, bin_edges_vec);
  
  cout << "c Marked ULC clauses" << endl;

  if (opts.orderencoderes) {

    // xlc_resolutions (ulcs);

    order_encode_resolutions (ulcs);

    // would only happend with general xlc resolution
    // if (unsat)
    //   return;

    // if (stats.order.unitsres) {
    //   // remove units from existing clauses (should not affect maps)
    //   // could mess up everything... maybe not...
    //   mark_satisfied_clauses_as_garbage();
    // }
  }

  cout << "c Performed resolutions" << endl;

  if (opts.orderencodealign) {

    // copy bin edges vec so we can delete them
    vector<set<int>> bin_edges_vec_copy;
    bin_edges_vec_copy.resize (2 * (vsize + 1));
    for (unsigned i = 0; i < 2 * (vsize + 1); i++) {
      bin_edges_vec_copy[i] = bin_edges_vec[i];
    }

    align_ulc_clauses (ulcs, bin_edges_vec_copy);

    cout << "c Aligned ULCs" << endl;
  }

  // Order Encoding (4 passes over clauses)
  //  1st, mark all clauses with variables having single occurence
  //       create noccs
  //  2nd, replace clause and resp. variables with order encoding
  //       check if clause meets criteria for order encoding
  //  3rd, use stamp to get rid of tautologies
  //  4th, replace literals in clauses with new literals
 
  // init_noccs ();

  // // Get occurence counts
  // LOG ("Get occurence counts");
  // for (auto c : clauses) {
  //   if (c->garbage) continue;
  //   for (auto lit : *c)
  //     noccs (lit)++;
  // }

  // get noccs (no need for occ lists!)
  vector<Clause*> garb_stack, proof_only;
  vector<int> pos_rep, neg_rep;
  for (int i = 0; i <= max_var; i++) {
    pos_rep.push_back (0);
    neg_rep.push_back (0);
    stamp[abs(i)] = 0;
  }

  int my_old_max_var = max_var;

  unordered_map<int,int> cliques;
  if (opts.orderencodereadclique)
    read_cliques (cliques);

  auto rng = std::default_random_engine {};

  vector<Clause*> pairwise_deletions;

  vector<int> var_in_hybrid;
  var_in_hybrid.resize ((vsize + 1), 0);

  LOG ("Find substitutions for order encoding");
  int cs = clauses.size (); // explicit size because we add clauses
  for (int c_idx = 0; c_idx < cs; c_idx++) {
    Clause *c = clauses[c_idx];
    if (c->garbage) continue; // possibly from resolutions
    // stats logging
    if (ulcs.find(c->id) != ulcs.end()) {

      if (opts.orderencodetype == 1) {
        if (ulcs[c->id] == 2 || ulcs[c->id] == 3 || ulcs[c->id] == 5)
          continue;
      }

      if (c->size <= opts.orderencodecut) {
        if (ulcs[c->id] == 1)
          stats.order.ulc_below++;
        else if (ulcs[c->id] == 2)
          stats.order.bin_below++;
        else if (ulcs[c->id] == 3)
          stats.order.hyb_below++;
        else if (ulcs[c->id] == 4)
          stats.order.both_below++;
        else if (ulcs[c->id] == 5)
          stats.order.hybf_below++;
          continue;
      }

      bool skipHybrid = false;
      if (ulcs[c->id] == 3) {
        for (auto lit: *c) {
          if (noccs (lit) == 1) {
            if (var_in_hybrid[abs(lit)] == 0) {
              var_in_hybrid[abs(lit)] = 1;
            } else {
              skipHybrid = true;
              break;
            }
          }
        }
      }
      if (skipHybrid) {
        continue;
      }

      if (ulcs[c->id] == 1)
        stats.order.nulc++;
      else if (ulcs[c->id] == 2)
        stats.order.nbin++;
      else if (ulcs[c->id] == 3)
        stats.order.nhyb++;
      else if (ulcs[c->id] == 4)
        stats.order.nboth++;
      else if (ulcs[c->id] == 5)
        stats.order.nhybf++;

      if (c->size > stats.order.max_size_reenc_clause)
        stats.order.max_size_reenc_clause = c->size;

      time += 1;
      LOG (c, "Can order clause");

      stats.order.foundorder++;
      vector<int> literals;
      vector<int> bin_literals;
      for (auto lit: *c) 
        literals.push_back (lit);

      

      bool is_bin = true;
      if (ulcs[c->id] == 1) {
        is_bin = false;
      }
      vector<Clause*> linear_encoding;
      if (!is_bin) {
        add_linear_encoding (literals, linear_encoding, c);
      } else if (ulcs[c->id] == 3) {
        // hybrid encoding
        add_pairwise_encoding (literals, bin_edges_vec, bin_edges_vec_cls, c);
      }
      
      // Becomes unsorted due to propagation in parsing
      // and/or preprocessing (e.g., lucky search)
      if (opts.orderencodeshuff) {
        shuffle (literals.begin(), literals.end(), rng);
      }
      else if (!opts.orderencodealign) { 
        // sort by absolute value
        sort (literals.begin(), literals.end(), [](int a, int b) {
          return abs(a) < abs(b);
        });
      }

      cout << "c order encoding: ";
      for (auto lit: literals) {
        cout << lit << " ";
      } cout << endl;

      int size = c->size;

      vector<vector<int>> temp_clauses;

      int aux = max_var + 1;
      external->internalize (aux); // all we need for new var?
      stats.order.new_vars++;

      // Define s0
      proof_only.push_back (add_proof_clause (vector<int>{-aux, literals[0]} , c, true));
      proof_only.push_back (add_proof_clause (vector<int>{aux, -literals[0]} , c, true));

      // perform variable elimination only if occ == 1
      if (opts.orderencodetype == 1 ||  noccs (literals[0]) == 1) {
        neg_rep[abs(literals[0])] = -aux;
        stamp[abs(literals[0])] = time;
        stats.order.nvars_ulc++;
      } else {
        bin_literals.push_back (literals[0]);
        stats.order.nvars_bin++;
      }

      for (int i = 1; i < size - 1; i++) {
        int v = abs(literals[i]);
        external->internalize (aux + 1);
        stats.order.new_vars++;

        // define si
        proof_only.push_back (add_proof_clause (vector<int>{-aux - 1, aux, literals[i]} , c, true));
        proof_only.push_back (add_proof_clause (vector<int>{aux + 1, -aux} , c, true));
        proof_only.push_back (add_proof_clause (vector<int>{aux + 1, -literals[i]} , c, true));
        proof_only.push_back (add_proof_clause (vector<int>{-aux, -literals[i]} , c, true));

        if (opts.orderencodetype == 1 || noccs (literals[i]) == 1) {
          pos_rep[v] = aux;
          neg_rep[v] = -aux - 1;
          stamp[abs(literals[i])] = time;
          stats.order.nvars_ulc++;
        } else {
          bin_literals.push_back (literals[i]);
          stats.order.nvars_bin++;
        }
        
        aux++;
      }

      // implicitly define sn (really with sn-1)
      proof_only.push_back (add_proof_clause (vector<int>{aux, literals[size-1]} , c, true));
      proof_only.push_back (add_proof_clause (vector<int>{-aux, -literals[size-1]} , c, true));
      
      if (opts.orderencodetype == 1 || noccs (literals[size-1]) == 1) {
        pos_rep[abs(literals[size-1])] = aux;
        stamp[abs(literals[size-1])] = time;
        stats.order.nvars_ulc++;
      } else {
        bin_literals.push_back (literals[size-1]);
        stats.order.nvars_bin++;
      }

      // temporarily mark garbage
      c->garbage = true;
      garb_stack.push_back(c);
      LOG (c, "Marked for deletion");

      // remove all bin_edges
      // if both ulc these will be removed as tautologies
      // but the stamp will not work for multiple occurrences of a literal
      if (is_bin) {
        pairwise_deletions.push_back(c); // must withhold deletions because some binary clauses may be used for multiple XLCs...
      } else {
        delete_linear_encoding (linear_encoding, my_old_max_var);
      }
    }
  }

  for (auto c : pairwise_deletions) {
    vector<int> literals;
    for (auto lit: *c) 
      literals.push_back (lit);
    delete_pairwise_encoding (literals, bin_edges_vec_cls);
  }

  cout << "c Finished defining order encodings" << endl;

  // stamp to delete clauses that will be tautologies (have repeated literals from reencoded clause)
  // not pefect, could be multiple stamps per clause...

  // Rup always by propagation of signals in order encoding
  // Prevents tautological clauses in case of exactly one constraints
  // only if there are not multiple containing same literal (update stamp)
  LOG ("Remove tautologies");
  for (auto c : clauses) {
    if (c->garbage) continue;
    bool del = false;
    int t = 0;
    for (auto lit : *c) {
      int s = stamp[abs(lit)]; 
      if (s) {
        if (s == t) {
          del = true;
          break;
        }
        if (!t)
          t = s;
      }
    }
    if (del) {
      LOG (c, "Tautology");
      stats.order.tautorder++;
      mark_garbage (c);
    }
  }

  // Reecode the clauses containing substituted literals.
  if (opts.orderencodeelim && (opts.orderencodetype == 1 || opts.orderencodetype == 3) && stats.order.nvars_ulc > 0) {
    // fill in for auxiliary variables added
    for (int i = pos_rep.size () ; i <= max_var; i++) {
      pos_rep.push_back (0);
      neg_rep.push_back (0);
    }
    stamp.resize ((vsize + 1) * 2, 0);

    LOG ("Reencode clauses");
    int s = clauses.size ();
    for (int i = 0; i < s; i++) {
      Clause *c = clauses[i];
      LOG (c, "Checking c for reencoding"); 
      if (c->garbage) continue;
      bool reencoding = false;
      for (auto lit : *c) {
        if (pos_rep[abs(lit)] != 0 || neg_rep[abs(lit)] != 0) {
          reencoding = true; // contains literal that was replaced
          break;
        }
      }

      if (reencoding) {

        if (c->size == 2)
          stats.order.nbin_subs++;

        time += 1;
        bool is_taut = false;
        LOG (c, "Reencode clause");
        for (auto lit : *c) {
          if (pos_rep[abs(lit)] != 0) {
            int sl = pos_rep[abs(lit)];
            if (stamp[vlit(-sl)] == time) {
              is_taut = true;
              break;
            }
            clause.push_back (sl);
            stamp[vlit(sl)] = time;
          } 
          if (neg_rep[abs(lit)] != 0) {
            int sl = neg_rep[abs(lit)];
            if (stamp[vlit(-sl)] == time) {
              is_taut = true;
              break;
            }
            clause.push_back (sl);
            stamp[vlit(sl)] = time;
          } 
          if (pos_rep[abs(lit)] == 0 && neg_rep[abs(lit)] == 0) {
            clause.push_back (lit);
          }
        }

        if (!is_taut) {
          stats.order.suborder++;
          new_clause_as (c);
        }
        else 
          stats.order.tautorder++;

        // Delete old clause. (Should be RUP)
        mark_garbage (c);

        clause.clear ();
      }
    }
    cout << "c finished substitutions" << endl;
  }

  
  vector<vector<int>> garb_stack_lits;

  // These are the ulc clauses to be deleted
  for (auto c : garb_stack) {
    c->garbage = false;
    vector<int> lits;
    for (auto lit : *c) {
      lits.push_back (lit);
    }
    garb_stack_lits.push_back(lits); // save to mark for elimination
    mark_garbage (c);
  }

  // delete proof only clauses
  bool keepAux = opts.orderencodekeep; // keep the clauses O1->O2 ... O(n-2)->O(n-1)

  // variable elimination removes definition clauses
  LOG ("Delete proof only");
  for (auto c : proof_only) { // unmark definition clauses
    c->garbage = false;
  }

  if (opts.orderencodeelim && (opts.orderencodetype == 1 || opts.orderencodetype == 3)) {
    for (auto c : proof_only) {      
      // clauses of the form si -> si+1
      if (c->size == 2 ) {
        if (abs(c->literals[0]) > my_old_max_var && abs(c->literals[1]) > my_old_max_var) {
          if (!keepAux)
            mark_garbage (c);
          continue;
        }
      }

      // definition clauses
      int pivot = c->literals[1];
      if (c->size > 2)
        pivot = c->literals[2];

      if (abs(pivot) > my_old_max_var) {
        cout << "ERROR: pivot in definition clause uses auxiliary variable" << endl;
        exit (1);
      }

      // bin case, no elimination
      if (pos_rep[abs(pivot)] == 0 && neg_rep[abs(pivot)] == 0)
        continue;

      // eliminate definition clauses
      if (c->size > 2)
        external->push_clause_on_extension_stack (c, pivot);
      else 
        external->push_binary_clause_on_extension_stack (c->id, pivot, c->literals[0]);
      
      mark_garbage (c);
    } 
  }

  cout << "c finished deletions" << endl;

  if (opts.orderencodeelim && (opts.orderencodetype == 1 || opts.orderencodetype == 3)) {
    // mark eliminated variables
    for (auto c : garb_stack_lits) {
      for (auto lit : c) {
        if (flags (abs(lit)).eliminated ())
          continue;
        if (pos_rep[abs(lit)] != 0 || neg_rep[abs(lit)] != 0) {
          stats.order.nvars_elim++;
          mark_eliminated (abs(lit));
        }
      }
    }
  }

  reset_noccs ();

  // alignment sorts literals inside a clause requiring watches to be reset
  if (stats.order.nres > 0 || stats.order.foundorder > 0 || opts.orderencodealign) {
    // get rid of binary garbage
    clear_watches ();
    connect_watches ();
    // reduce (); // process units from resolutions?
    // binary clauses not deleted until next call to reduce, may not appear in proof if 
    // problem is solved before then
  }

  LOG ("STOP order encoding");

  STOP (order);

  // if (stats.order.nulc + stats.order.nbin + stats.order.nhyb + stats.order.nboth + stats.order.nhybf == 0) {
  //   // internal->stats.print (internal);
  //   cout << "c Nothing Found Exiting" << endl;
  //   internal->stats.print (internal);
  //   exit (1);
  // }

  // if (opts.orderencodetype == 3 && stats.order.nbin + stats.order.nhyb + stats.order.nhybf == 0) {
  //   cout << "c Nothing NonULC Found Exiting" << endl;
  //   internal->stats.print (internal);
  //   exit (1);
  // }

  if (opts.orderencodeexit) {
    cout << "c Exiting after order encoding" << endl;
    // first print statistics
    internal->stats.print (internal);
    exit (1);
  }

}

bool order_encoding () {
  // when to perform order encoding (always)
  return true;
}

} // namespace CaDiCaL
