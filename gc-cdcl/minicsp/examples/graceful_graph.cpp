/*************************************************************************
minicsp

Copyright 2010 Nina Narodytska

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
#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"
#include "minicsp/core/cmdline.hpp"
#include "minicsp/core/utils.hpp"

using namespace std;
using namespace minicsp;

int main(int argc, char *argv[])
{
  cmdline::arglist args(argv+1, argv+argc);
  Solver s;
  cmdline::parse_solver_options(s, args);
  bool stat = cmdline::has_option(args, "--stat");
  setup_signal_handlers(&s);

  if( args.empty() ) {
    cerr << "usage " << argv[0] << " [options] <size>\n";
    return 1;
  }

  size_t g = atoi(args.back().c_str());
  int n = 2 * g + 1;
  int q = 2 * g + 2 * g;

  vector<cspvar> x = s.newCSPVarArray(n, 0, q);
  vector<cspvar> diff = s.newCSPVarArray(q, -q, q);
  vector<cspvar> abs_diff = s.newCSPVarArray(q, 1, q);


  post_alldiff(s, x);
  post_alldiff(s, abs_diff);


  // the first circle
  int k = 0;
  for(size_t i = 0; i < g; ++i) {
      vector<cspvar> v(3);
      vector<int> c(3);
      v[0] = x[i]; c[0] = -1;
      v[1] = x[(i + 1) % g]; c[1] = 1;
      v[2] = diff[k]; c[2] = -1;
      post_lin_eq(s, v, c, 0);
      k++;
    }
  // the second circle
  for(size_t i = 0; i < g; ++i) {
      vector<cspvar> v(3);
      vector<int> c(3);
      v[0] = x[g + i]; c[0] = -1;
      v[1] = x[g + (i + 1) % g]; c[1] = 1;
      v[2] = diff[k]; c[2] = -1;
      post_lin_eq(s, v, c, 0);
      k++;
    }
  // connections
  for(size_t i = 0; i < 2 * g; ++i) {
      vector<cspvar> v(3);
      vector<int> c(3);
      v[0] = x[i]; c[0] = -1;
      v[1] = x[2 * g]; c[1] = 1;
      v[2] = diff[k]; c[2] = -1;
      post_lin_eq(s, v, c, 0);
      k++;
  }
  cout << " k " << k - 1 << " q " << q << endl;
  for(int i = 0; i < q; ++i) {
          post_abs(s, diff[i], abs_diff[i], 0);
  }

 // circle 1: x[0] < others
  for(size_t i = 1; i < g; ++i)
    post_less(s, x[0], x[i], 0);
// +
  post_less(s, x[1], x[g - 1], 0);

// circle 2: x[g] < others
  for(size_t i = 1; i < g; ++i)
    post_less(s, x[g], x[i + g], 0);
// +
  post_less(s, x[1 + g], x[g + g - 1], 0);

// between circles
    post_less(s, x[0], x[g], 0);

// inner symmetry
    x[2 * g].assign(s, 0, NO_REASON);

// circle 1: x[0]  = 1 x[2] = 3 ...
    for(size_t i = 0; i < g - 3; i = i + 2)
     x[i].assign(s, i + 1, NO_REASON);

// circle 2: x[0]  = 2 x[2] = 4 ...
   for(size_t i = 0; i < g - 3; i = i + 2)
        x[i + g].assign(s, i + 2, NO_REASON);

  bool sol = false;
  sol = s.solve();

  cout << s.conflicts << " conflicts\n";
  if( !sol ) {
    cout << "unsat\n";
    return 0;
  }
  cout << "the first circle " << endl;
  for(size_t i = 0; i < g; ++i) {
          cout << " n [" <<  i << "] " << s.cspModelValue(x[i]) << "; ";
  }
  cout << endl;
  cout << "the second circle " << endl;
  for(size_t i = 0; i < g; ++i) {
          cout << " n [" <<  g + i << "] " << s.cspModelValue(x[g + i]) << "; ";
    }
  cout << endl;

  k = 0;
  cout << "diff: the first circle " << endl;
  for(size_t i = 0; i < g; ++i) {
          cout << " n [" <<  i << "," <<  (i + 1) % g  <<  "] = |" << s.cspModelValue(x[i]) - s.cspModelValue(x[(i + 1) % g]) << "| =  "<< s.cspModelValue(abs_diff[k]) << "; ";                  k++;
  }
  cout << endl;
  cout << "diff: the second circle " << endl;
    for(size_t i = 0; i < g; ++i) {
                cout << " n [" <<  i + g << "," << g + (i + 1) % g  <<  "] = |" << s.cspModelValue(x[i + g]) - s.cspModelValue(x[(i + 1) % g + g]) << "| =  " << s.cspModelValue(abs_diff[k]) << "; ";
                  k++;
    }
    cout << endl;


    cout << "diff: connections " << endl;
    for(size_t i = 0; i < 2*g; ++i) {
                cout << " n [" <<  2*g << "," <<  i <<  "] = |" << s.cspModelValue(x[2*g]) - s.cspModelValue(x[i]) << "| =  " << s.cspModelValue(abs_diff[k]) << "; ";
                  k++;
    }
    cout << endl;

  if( stat )
    printStats(s);

  return 0;
}
