/*************************************************************************
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

*************************************************************************/

/* SEND+MORE=MONEY */

#include <iostream>
#include <vector>
#include "minicsp/core/solver.hpp"
#include "minicsp/core/cons.hpp"

using namespace std;
using namespace minicsp;

int main()
{
  Solver s;
  vector<cspvar> x = s.newCSPVarArray(8, 0, 9);
  cspvar S = x[0];
  cspvar E = x[1];
  cspvar N = x[2];
  cspvar D = x[3];
  cspvar M = x[4];
  cspvar O = x[5];
  cspvar R = x[6];
  cspvar Y = x[7];

  vector<cspvar> v;
  vector<int> c;
  v.push_back(S); c.push_back(1000);
  v.push_back(E); c.push_back(100);
  v.push_back(N); c.push_back(10);
  v.push_back(D); c.push_back(1);
  v.push_back(M); c.push_back(1000);
  v.push_back(O); c.push_back(100);
  v.push_back(R); c.push_back(10);
  v.push_back(E); c.push_back(1);
  v.push_back(M); c.push_back(-10000);
  v.push_back(O); c.push_back(-1000);
  v.push_back(N); c.push_back(-100);
  v.push_back(E); c.push_back(-10);
  v.push_back(Y); c.push_back(-1);

  vector<int> c1(c);
  for(size_t i = 0; i != c1.size(); ++i) c1[i] = -c1[i];

  post_lin_leq(s, v, c, 0);
  post_lin_leq(s, v, c1, 0);

  s.solve();
}
