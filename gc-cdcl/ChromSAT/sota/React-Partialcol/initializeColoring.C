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

#include "initializeColoring.h"


void initializeColoring(Graph & g, int * c, int k, Random & r) {

    // A simple greedy algorithm that leaves the assigned color
    // if possible, gives another legal color or
    // assigns color 0 if nothing
    // else is available.

    // First produce a random permutation for the vertex order
    int perm[g.n];
    for (int i = 0; i < g.n; i++) {
        perm[i] = i;
    }
    for (int i = 0; i < g.n; i++) {
        int p = r.getInt(0, g.n - 1);
        int h = perm[i];
        perm[i] = perm[p];
        perm[p] = h;
    }

    int taken[k + 1];

    // Insure all colors are in the range [0, ... ,k]
    for (int i = 0; i < g.n; i++) {
        if (c[i] < 0 || c[i] > k)
            c[i] = 0;
    }

    // Go through all nodes
    for (int ii = 0; ii < g.n; ii++) {
        int i = perm[ii];
        // Build a list of used colors in the nodes neighborhood
        for (int j = 0; j <= k; j++) {
            taken[j] = 0;
        }
        // for (int j = 0; j < g.n; j++) {
        //     if (i != j && g[i][j]) {
        //         taken[c[j]]++;
        //     }
        // }
        for (std::vector<int>::iterator j = g.neighbor[i].begin();
             j != g.neighbor[i].end(); j++) {
            taken[c[*j]]++;
        }
        // if the currently assigned color is legal and not 0, leave it
        // otherwise find a new legal color, and if not possible
        // set it to zero.
        if (c[i] == 0 || taken[c[i]] > 0) {
            int color = 0;
            for (int j = 0; j <= k; j++) {
                if (taken[j] == 0) {
                    color = j;
                    break;
                }
            }
            c[i] = color;
        }
    }
}







void initializeColoringForTabu(Graph & g, int * c, int k, Random & r) {

    // A simple greedy algorithm that leaves the assigned color
    // if possible, gives another legal color or
    // assigns a random color if nothing
    // else is available.

    // First produce a random permutation for the vertex order
    int perm[g.n];
    for (int i = 0; i < g.n; i++) {
        perm[i] = i;
    }
    for (int i = 0; i < g.n; i++) {
        int p = r.getInt(0, g.n - 1);
        int h = perm[i];
        perm[i] = perm[p];
        perm[p] = h;
    }

    int taken[k + 1];

    // Insure all colors are in the range [1, ... ,k]
    for (int i = 0; i < g.n; i++) {
        if (c[i] < 1 || c[i] > k)
            c[i] = 1;
    }

    // Go through all nodes
    for (int ii = 0; ii < g.n; ii++) {
        int i = perm[ii];
        // Build a list of used colors in the nodes neighborhood
        for (int j = 1; j <= k; j++) {
            taken[j] = 0;
        }
        // for (int j = 0; j < g.n; j++) {
        //     if (i != j && g[i][j]) {
        //         taken[c[j]]++;
        //     }
        // }
        for (std::vector<int>::iterator j = g.neighbor[i].begin();
             j != g.neighbor[i].end(); j++) {
            taken[c[*j]]++;
        }
        // if the currently assigned color is legal, leave it
        // otherwise find a new legal color, and if not possible
        // set it to a random color.
        if (taken[c[i]] > 0) {
            int color = r.getInt(1, k);
            for (int j = 1; j <= k; j++) {
                if (taken[j] == 0) {
                    color = j;
                    break;
                }
            }
            c[i] = color;
        }
    }
}






