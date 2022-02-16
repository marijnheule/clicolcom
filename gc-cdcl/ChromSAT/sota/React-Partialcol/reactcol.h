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

#ifndef REACTCOL_INCLUDED
#define REACTCOL_INCLUDED

#include "Graph.h"
#include "Random.h"

int reactcol(Graph & g, int * c, int k, Random & r,
	     int maxIterations, int tenure, int verbose, 
	     int frequency, int increment, int maxTime);


#endif
