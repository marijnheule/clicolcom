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

#include "manipulateArrays.h"

#include <iostream>
#include "termcolors.h"

void initializeArrays(int**& nodesByColor, int**& conflicts, int**& tabuStatus,
    int**& neighbors, int*& nbcPosition, Graph& g, int* c, int k)
{
    int n = g.n;

    // Allocate and initialize (k+1)x(n+1) array for nodesByColor and conflicts
    nodesByColor = new int*[k + 1];
    conflicts = new int*[k + 1];
    for (int i = 0; i <= k; i++) {
        nodesByColor[i] = new int[n + 1];
        nodesByColor[i][0] = 0;
        conflicts[i] = new int[n + 1];
        for (int j = 0; j <= n; j++) {
            conflicts[i][j] = 0;
        }
    }

    // Allocate the neighbors array
    neighbors = new int*[n];
    for (int i = 0; i < n; i++) {
        neighbors[i] = new int[n + 1];
        neighbors[i][0] = 0;
    }

    // Allocate the tabuStatus array
    tabuStatus = new int*[n];
    for (int i = 0; i < n; i++) {
        tabuStatus[i] = new int[k + 1];
        for (int j = 0; j <= k; j++) {
            tabuStatus[i][j] = 0;
        }
    }

    // Allocate the nbcPositions array
    nbcPosition = new int[n];

    // Initialize the nodesByColor and nbcPosition array
    for (int i = 0; i < n; i++) {
        // C is cool ;-)
        nodesByColor[c[i]][(nbcPosition[i] = ++nodesByColor[c[i]][0])] = i;
    }

    // Initialize the conflicts and neighbors array
    for (int i = 0; i < n; i++) {
        // for (int j = 0; j < n; j++) {
        //     if (g[i][j] && i != j) {
        //         conflicts[c[j]][i]++;
        //         neighbors[i][++neighbors[i][0]] = j;
        //     }
        // }
      for (std::vector<int>::iterator j = g.neighbor[i].begin(); j != g.neighbor[i].end(); j++) {
              conflicts[c[*j]][i]++;
              neighbors[i][++neighbors[i][0]] = *j;
      }
    }
}

void moveNodeToColor(int bestNode, int bestColor, Graph& g, int* c,
    int** nodesByColor, int** conflicts, int* nbcPosition, int** neighbors,
    int** tabuStatus, int totalIterations, int tabuTenure)
{

    // move bestNodes to bestColor
    c[bestNode] = bestColor;
    // Replace bestNode by the last node in the nodesByColor[0] array and
    // shorten it
    nodesByColor[0][nbcPosition[bestNode]]
        = nodesByColor[0][nodesByColor[0][0]--];
    // Update the nbcPosition array the node that has taken the place of best
    // node
    nbcPosition[nodesByColor[0][nbcPosition[bestNode]]] = nbcPosition[bestNode];
    // Insert bestNode into the nodesByColor[bestColor] array, increase its
    // lenght
    // and update the nbcPosition for bestNode
    nodesByColor[bestColor]
                [(nbcPosition[bestNode] = ++nodesByColor[bestColor][0])]
        = bestNode;

    //  cout << "moved node " << NODE(bestNode) << " to color " <<
    //  COLOR(bestColor) << endl;

    // Update the conflicts array and remove conflicting nodes
    for (int j = 1; j <= neighbors[bestNode][0]; j++) {
        int i = neighbors[bestNode][j];

        // Do not move neighbors to bestColor for a couple of iterations in
        // order to
        // avoid bestNode from dropping back out too soon
        tabuStatus[i][bestColor] = totalIterations + tabuTenure;

        // Increase the conflicts for bestColor
        conflicts[bestColor][i]++;
        // Check for conflict created by moving bestNode to bestColor


        if (c[i] == bestColor) {
            // If so, remove the node and put it back to the uncolored nodes
            nodesByColor[bestColor][nbcPosition[i]]
                = nodesByColor[bestColor][nodesByColor[bestColor][0]--];
            nbcPosition[nodesByColor[bestColor][nbcPosition[i]]]
                = nbcPosition[i];
            nodesByColor[0][(nbcPosition[i] = ++nodesByColor[0][0])] = i;
            c[i] = 0;
            // Reduce the conflicts of all neighbors.
            for (int k = 1; k <= neighbors[i][0]; k++) {
                conflicts[bestColor][neighbors[i][k]]--;
            }
            //      cout << "removed node " << NODE(i) << " to the OUT set\n";
        }
    }
}

void moveNodeToColorForTabu(int bestNode, int bestColor, Graph& g, int* c,
    int** nodesByColor, int** conflicts, int* nbcPosition, int** neighbors,
    int* nodesInConflict, int* confPosition, int** tabuStatus,
    int totalIterations, int tabuTenure)
{

    int oldColor = c[bestNode];
    // move bestNodes to bestColor
    c[bestNode] = bestColor;

    // If bestNode is not a conflict node anymore, remove it from the list
    if (conflicts[oldColor][bestNode] && !(conflicts[bestColor][bestNode])) {
        confPosition[nodesInConflict[nodesInConflict[0]]]
            = confPosition[bestNode];
        nodesInConflict[confPosition[bestNode]]
            = nodesInConflict[nodesInConflict[0]--];
    } else // If bestNode becomes a conflict node, add it to the list
        if (!(conflicts[oldColor][bestNode])
            && conflicts[bestColor][bestNode]) {
        nodesInConflict[(confPosition[bestNode] = ++nodesInConflict[0])]
            = bestNode;
    }

    // Update the conflicts of the neighbors.
    for (int i = 1; i <= neighbors[bestNode][0]; i++) {
        int nb = neighbors[bestNode][i];
        // Decrease the number of conflicts in the old color
        if ((--conflicts[oldColor][nb]) == 0 && c[nb] == oldColor) {
            // Remove nb from the list of conflicting nodes if there are 0
            // conflicts in
            // its own color
            confPosition[nodesInConflict[nodesInConflict[0]]]
                = confPosition[nb];
            nodesInConflict[confPosition[nb]]
                = nodesInConflict[nodesInConflict[0]--];
        }
        // Increase the number of conflicts in the new color
        if ((++conflicts[bestColor][nb]) == 1 && c[nb] == bestColor) {
            // Add nb from the list conflicting nodes if there is a new conflict
            // in
            // its own color
            nodesInConflict[(confPosition[nb] = ++nodesInConflict[0])] = nb;
        }
    }
    // Set the tabu status
    tabuStatus[bestNode][oldColor] = totalIterations + tabuTenure;
}

void freeArrays(int**& nodesByColor, int**& conflicts, int**& tabuStatus,
    int**& neighbors, int*& nbcPosition, int k, int n)
{

    for (int i = 0; i <= k; i++) {
        delete[] nodesByColor[i];
        delete[] conflicts[i];
    }
    for (int i = 0; i < n; i++) {
        delete[] neighbors[i];
        delete[] tabuStatus[i];
    }
    delete[] nodesByColor;
    delete[] conflicts;
    delete[] tabuStatus;
    delete[] nbcPosition;
    delete[] neighbors;
}

