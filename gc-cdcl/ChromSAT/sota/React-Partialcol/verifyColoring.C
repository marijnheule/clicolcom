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

#include "verifyColoring.h"

#include <iostream>
#include "termcolors.h" 

using namespace std;

int verifyColoring(Graph & g, int * coloring, int k) {

	for (int i=0; i<g.n; i++) {
		if (coloring[i] < 1 || coloring[i] > k) {

			cout << ERROR("color out of range") << " for node " 
				<< NODE(i) << " (" << COLOR(coloring[i]) << ")\n";
			return 0;
		}
		for(std::vector<int>::iterator j = g.neighbor[i].begin(); j != g.neighbor[i].end(); ++j) {
			if (coloring[i] == coloring[*j]) {
				cout << ERROR("Conflict!") << " between node " << NODE(i) << " and " << NODE(*j) << endl;
				return 0;
			}
		}
		// for (int j=i+1; j<g.n; j++) {
		// 	if (g[i][j]) {
		// 		if (coloring[i] == coloring[j]) {
		// 			cout << ERROR("Conflict!") << " between node " << NODE(i) << " and " << NODE(j) << endl;
		// 			return 0;
		// 		}
		// 	}
		// }
	}
	return 1;
}
