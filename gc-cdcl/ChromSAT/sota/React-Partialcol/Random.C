/******************************************************************************/
//
//  ReactPartialCol, PartialCol, ReactTabucol and Tabucol graph coloring
//  heuristics. Reference code for the paper
//  "A Reactive Tabu Search Using Partial Solutions for the
//  Graph Coloring Problem" by Ivo Bloechliger and Nicolas Zuffery.
//
//  Copyright (C) 2003 - Ecole Poyltechnique Federale de Lausanne - EPFL, 
//  Laboratory of Operations Research South Est-ROSE, 1015 Lausanne, Switzerland
//  Written by Ivo Bloechliger, Ivo.Bloechliger@epfl.ch
//  http://rose.epfl.ch/~bloechli/coloring/
//
/******************************************************************************/
//
// This program is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation. In paticular, this program is 
// distributed WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The EPFL shall in no 
// case be liable for any damage of any kind in connection with the use of this
// program.  See the GNU General Public License for more details 
// (http://www.gnu.org/copyleft/gpl.html#SEC1).
//
/******************************************************************************/


#include "Random.h"

Random :: Random(int aSeed) {
  this->seed=aSeed;
}

int Random :: getInt(int a1, int a2) {
  return (int) (a1 + this->getDouble() * (a2-a1+1));

}

double Random :: getDouble() {
  long m, a, b, c;
  long k1;
  double tempo;
  
  m = 2147483647L;
  a = 16807;
  b = 127773L;
  c = 2836;
  
  k1 = this->seed / b;
  this->seed = a * (this->seed % b) - k1 * c;
  if (this->seed < 0) this->seed += m;
  tempo = ((double) this->seed) / ((double) m);
  if (tempo == 0) tempo = 0.5;
  return tempo;
}

void Random::permutation(int * perm, int n) {
  for (int i=n-1; i>=0; i--) {
    int pos = getInt(0,i);
    int v=perm[pos];
    perm[pos] = perm[i];
    perm[i] = v;
  }
}

