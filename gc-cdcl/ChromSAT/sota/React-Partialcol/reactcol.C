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

#include "reactcol.h"

#include "initializeColoring.h" 
#include "manipulateArrays.h"

#include <iostream>
#include "termcolors.h"

#include "checkCPUtime.h"


using namespace std;

int reactcol(Graph& g, int* c, int k, Random& r, int maxIterations,
    int staticTenure, int verbose, int freq, int inc, int maxTime)
{

    int** nodesByColor; // Arrays of nodes for each color
    int* nbcPosition; // Position of each node in the above array
    int** conflicts; // Number of conflicts for each color and node
    int** tabuStatus; // Tabu status for each node and color
    int** neighbors; // List of neighbors for each node

    //  int pairs[][3] = {{10000,10,10}};

    int pairs[][3]
        = {{10000, 10, 5}, {10000, 15, 3}, {10000, 5, 10}, {5000, 15, 10},
            {5000, 10, 15}, {5000, 5, 20}, {1000, 15, 30}, {1000, 10, 50},
            {1000, 5, 100}, {500, 5, 100}, {500, 10, 150}, {500, 15, 200}};

    int numPairs = sizeof(pairs) / sizeof(int) / 3;

    int pairCycles = 0;
    int frequency = pairs[0][0];
    int increment = pairs[0][1];
    int nextPair = pairs[0][2];

    if (freq) { // frenquency and increment are set manually
        frequency = freq;
        increment = inc;
    }

    int totalIterations = 0;
    int currentIterations = 0;
    int incVerbose = 0;
    if (verbose == 1)
        incVerbose = 10000;
    if (verbose == 2)
        incVerbose = 100;
    if (verbose > 2)
        incVerbose = 1;
    int nextVerbose = incVerbose;
    int result = -1;

    int tabuTenure = increment;
    int randomTenure = 1;
    if (staticTenure != 0) {
        tabuTenure = staticTenure;
        randomTenure = 0;
    }

    initializeColoring(g, c, k, r);
    if (verbose > 1)
        cout << "Initialized the coloring\n";

    initializeArrays(
        nodesByColor, conflicts, tabuStatus, neighbors, nbcPosition, g, c, k);
    if (verbose > 1)
        cout << "Initialized the arrays. |Outnodes| = " << nodesByColor[0][0]
             << endl;

    int bestSolutionValue = nodesByColor[0][0]; // Number of out nodes

    // Just in case we already have an admissible k-coloring
    if (bestSolutionValue == 0) {
        return 0;
    }

    int minSolutionValue = g.n;
    int maxSolutionValue = 0;

    int numRandDec = 0;

    while (currentIterations < maxIterations && checkCPUtime() < maxTime) {

        currentIterations++;
        totalIterations++;

        int bestNode = -1, bestColor = -1, bestValue = g.n;
        int numBest = 0;

        // Try for every uncolored outNode

        for (int iOutNode = 1; iOutNode <= nodesByColor[0][0]; iOutNode++) {
            int outNode = nodesByColor[0][iOutNode];
            // to move it to every color

            for (int color = 1; color <= k; color++) {
                if (conflicts[color][outNode] <= bestValue) {
                    if (conflicts[color][outNode] < bestValue) {
                        numBest = 0;
                    }

                    // Only consider the move if it is not tabu or leads
                    // to a new very best solution seen globally.

                    if (tabuStatus[outNode][color] < totalIterations
                        || (conflicts[color][outNode] == 0
                               && nodesByColor[0][0] == bestSolutionValue)) {

                        // Select the nth move with probability 1/n

                        if (r.getInt(0, numBest) == 0) {
                            bestNode = outNode;
                            bestColor = color;
                            bestValue = conflicts[color][outNode];
                        }
                        numBest++; // Count the number of considered moves
                    }
                }
            }
        }
        // If no non tabu moves have been found, take any random move
        if (bestNode == -1) {
            ++numRandDec;
            bestNode = nodesByColor[0][r.getInt(1, nodesByColor[0][0])];
            bestColor = r.getInt(1, k);
            bestValue = conflicts[bestColor][bestNode];
        }

        int tTenure = tabuTenure;
        if (randomTenure == 1)
            tTenure = r.getInt(1, tTenure);
        // Now execute the move
        moveNodeToColor(bestNode, bestColor, g, c, nodesByColor, conflicts,
            nbcPosition, neighbors, tabuStatus, totalIterations, tTenure);

        if (verbose > 2) {
            cout << "Moved node " << NODE(bestNode) << " to color "
                 << COLOR(bestColor) << " with value " << NUM(bestValue)
                 << endl;
        }

        // Update the min and max objective function value
        if (nodesByColor[0][0] > maxSolutionValue)
            maxSolutionValue = nodesByColor[0][0];
        if (nodesByColor[0][0] < minSolutionValue)
            minSolutionValue = nodesByColor[0][0];

        int Delta = maxSolutionValue - minSolutionValue;

        if (staticTenure == 0) {

            if (currentIterations % frequency == 0) {
                // Adjust the tabuTenure every frequency iterations
                if (Delta < 2 || tabuTenure == 0) {
                    tabuTenure += increment;
                    if (pairCycles == nextPair) {
                        if (!freq) { // frequency and incrment are not set
                            // manually
                            int p = r.getInt(0, numPairs - 1);
                            frequency = pairs[p][0];
                            increment = pairs[p][1];
                            pairCycles = 0;
                            nextPair = pairs[p][2];
                        }
                        randomTenure = r.getInt(0, 1);
                    } else {
                        pairCycles++;
                    }
                } else if (tabuTenure) {
                    tabuTenure--;
                }

                minSolutionValue = g.n;
                maxSolutionValue = 0;
            }
        } else {
            tabuTenure = (int)(0.6 * nodesByColor[0][0]) + r.getInt(0, 9);
        }

        // Have we a new globally best solution?
        if (nodesByColor[0][0] < bestSolutionValue) {
            bestSolutionValue = nodesByColor[0][0];

            // If all nodes are colored we report success and stop iterating
            if (bestSolutionValue == 0) {
                result = 1;
                break;
            }
            // Otherwise reinitialize some values
            minSolutionValue = g.n;
            maxSolutionValue = 0;
            currentIterations = 0;
            pairCycles = 0;
            nextVerbose = totalIterations;
        }
        if (totalIterations == nextVerbose && incVerbose) {
            cout << NUMBER(totalIterations) << NUMBER(currentIterations)
                 << "   obj =" << NUM(nodesByColor[0][0])
                 << "   best =" << NUM(bestSolutionValue)
                 << "   tenure =" << NUM(tabuTenure)
                 << "   freq = " << NUMBER(frequency)
                 << "   inc = " << NUM(increment)
                 << "   rand = " << NUM(randomTenure)
                 << "   Delta = " << NUM(Delta) << endl;
            nextVerbose += incVerbose;
        }
    }

    freeArrays(
        nodesByColor, conflicts, tabuStatus, neighbors, nbcPosition, k, g.n);

    return result * totalIterations;
}




