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

#include "tabu.h"
#include "initializeColoring.h" 
#include "manipulateArrays.h"

#include <iostream>
#include "termcolors.h"

#include "checkCPUtime.h"

using namespace std;

int tabu(Graph & g, int * c, int k, Random & r,
	 int maxIterations, int staticTenure, int verbose, 
	 int freq, int inc, int maxTime) {

  int ** nodesByColor; // Arrays of nodes for each color
  int * nbcPosition;   // Position of each node in the above array
  int ** conflicts;   // Number of conflicts for each color and node
  int ** tabuStatus;  // Tabu status for each node and color
  int ** neighbors;   // List of neighbors for each node
  int nodesInConflict[g.n+1];
  int confPosition[g.n];

  int pairs[][3] = {{10000,10,5},
		    {10000,15,3},
		    {10000,5,10},
		    {5000,15,10},
		    {5000,10,15},
		    {5000,5,20},
		    {1000,15,30},
		    {1000,10,50},
		    {1000,5,100},
		    {500,5,100},
		    {500,10,150},
		    {500,15,200}};

  int numPairs = sizeof(pairs)/sizeof(int)/3;

  int pairCycles = 0;
  int frequency = pairs[0][0];
  int increment = pairs[0][1];
  int nextPair = pairs[0][2];

  if (freq) {
    frequency = freq;
    increment = inc;
  }


  int totalIterations = 0;
  int currentIterations = 0;
  int incVerbose = 0;
  int totalConflicts = 0;
  if (verbose == 1) incVerbose=10000;
  if (verbose == 2) incVerbose=100;
  if (verbose > 2) incVerbose=1;
  int nextVerbose = incVerbose;
  int result = -1;

  int tabuTenure = g.n/10;
  int randomTenure = 1;
  if (staticTenure != 0) {
    tabuTenure = staticTenure;
    randomTenure = 0;
  }

  initializeColoringForTabu(g, c, k, r);

  if (verbose>1) cout << "Initialized the coloring\n";
  
  initializeArrays(nodesByColor, conflicts, tabuStatus, 
		   neighbors, nbcPosition, g, c, k);
  // Count the number of conflicts and set up the list nodesInConflict
  // with the associated list confPosition
  nodesInConflict[0]=0;
  for (int i=0; i<g.n; i++) {
    if (conflicts[c[i]][i] > 0) {
      totalConflicts += conflicts[c[i]][i];
      nodesInConflict[ (confPosition[i]=++nodesInConflict[0]) ] = i;
    }   
  }
  totalConflicts /=2;

  if (verbose>1) cout << "Initialized the arrays. #Conflicts = " << totalConflicts << endl;
  
  int bestSolutionValue = totalConflicts; // Number of conflicts

  // Just in case we already have an admissible k-coloring
  if (bestSolutionValue == 0) {
    return 0;
  }

  int minSolutionValue = g.n;
  int maxSolutionValue = 0;

  while (currentIterations < maxIterations && checkCPUtime() < maxTime) {

    currentIterations++;
    totalIterations++;

    int nc = nodesInConflict[0];

    int bestNode=-1, bestColor=-1, bestValue=g.n*g.n;
    int numBest=0;

//      int cnt=0;
//      for (int i=0; i<g.n; i++) {
//        if (conflicts[c[i]][i]) {
//  	cnt++;
//  	if (nodesInConflict[confPosition[i]] != i) {
    //  	  cout << "Problem with node " << NODE(i) << endl;
//  	}
//        }
//      }
//      if (cnt != nodesInConflict[0]) {
//        cout << "nodesInConflicts count is wrong! real=" << cnt << " but nic = " << nodesInConflict[0] << endl;
//      }
    
    // Try for every node in conflict
    
    for (int iNode=1; iNode <= nodesInConflict[0]; iNode++) {
      int node = nodesInConflict[iNode];
      // to move it to every color

      for (int color=1; color<=k; color++) {
	if (color != c[node]) {
	  int newValue = totalConflicts + conflicts[color][node] - conflicts[c[node]][node];
//  	  cout  << totalConflicts << " + " <<  conflicts[color][node] << " - " << conflicts[c[node]][node] 
//  		<< " wich makes " << newValue << endl;
	  if (newValue <= bestValue && color != c[node]) {
	    if (newValue < bestValue) {
	      numBest=0;
	    }
//  	    cout << " ---- Found newValue " << newValue << " with node " 
//  		 << NODE(node) << " and color " << COLOR(color) << endl;
	    // Only consider the move if it is not tabu or leads
	    // to a new very best solution seen globally.
	    
	    if (tabuStatus[node][color] < totalIterations || 
		(newValue < bestSolutionValue)) {	  
	      
	      // Select the nth move with probability 1/n
	      
	      if (r.getInt(0,numBest)==0) {
		bestNode = node;
		bestColor = color;
		bestValue = newValue;
	      }
	      numBest++;  // Count the number of considered moves
	    }
	  }
	}
      }
    }
    // If no non tabu moves have been found, take any random move
    if (bestNode == -1) {
      bestNode = r.getInt(0, g.n-1);
      while ((bestColor = r.getInt(1,k)) != c[bestNode]);
      bestValue = totalConflicts + conflicts[bestColor][bestNode] - conflicts[c[bestNode]][bestNode];
    }
    
    // Now execute the move
    if (verbose>2) {
      cout << "Will move node " << NODE(bestNode) << " to color " << COLOR(bestColor) 
	   << " with value " << NUM(bestValue) 
	   << " oldconf = " << NUM(conflicts[c[bestNode]][bestNode])
	   << " newconf = " << NUM(conflicts[bestColor][bestNode])
	   << " totalConflicts = " << NUM(totalConflicts)
	   << endl;
    }
    int tTenure = tabuTenure;
    if (randomTenure == 1) tTenure = r.getInt(1,tTenure);
    // Now execute the move
    moveNodeToColorForTabu(bestNode, bestColor, g, c, 
			   nodesByColor, conflicts, 
			   nbcPosition, neighbors,
			   nodesInConflict, confPosition,
			   tabuStatus, totalIterations, tTenure);
    totalConflicts = bestValue;

    int max_min = 0;

    if (staticTenure == 0) {
      // Update the min and max objective function value
      if (totalConflicts > maxSolutionValue) maxSolutionValue = totalConflicts;
      if (totalConflicts < minSolutionValue) minSolutionValue = totalConflicts;
      
      max_min = maxSolutionValue - minSolutionValue;
     
      if (currentIterations % frequency == 0) {
	// Adjust the tabuTenure every frequency iterations
	if (maxSolutionValue - minSolutionValue < 4+totalConflicts/80 || tabuTenure==0) {
	  tabuTenure += increment;
	  if (pairCycles == nextPair) {
	    if (!freq) {// frequency and increment are not set manually
	      int p = r.getInt(0,numPairs-1);
	      frequency = pairs[p][0];
	      increment = pairs[p][1];
	      pairCycles = 0;
	      nextPair = pairs[p][2];
	    }
	    randomTenure = r.getInt(0,1);
	  }
	} else if (tabuTenure) {
	  tabuTenure--;
	}
	
	minSolutionValue = g.n*g.n;
	maxSolutionValue = 0;
	
	if (pairCycles == nextPair) {
	  if (!freq) { // frequency and increment are not set manually
	    int p = r.getInt(0,numPairs-1);
	    frequency = pairs[p][0];
	    increment = pairs[p][1];
	    pairCycles = 0;
	    nextPair = pairs[p][2];
	  }
	} else {
	  pairCycles++;
	}
	
      }
    } else {
      tabuTenure = (int)(0.6*nc) + r.getInt(0,9);
    }
    
    // Have we a new globally best solution?
    if (totalConflicts < bestSolutionValue) {
      bestSolutionValue = totalConflicts;
      
      // If all nodes are colored we report success and stop iterating
      if (bestSolutionValue == 0) {
	result = 1;
	break;
      }
      // Otherwise reinitialize some values
      minSolutionValue = g.n*g.n;
      maxSolutionValue = 0;
      currentIterations = 0;
      nextVerbose = totalIterations;
    }
    if (totalIterations == nextVerbose && incVerbose) {
      cout << NUMBER(totalIterations) << NUMBER(currentIterations)
	   << "   obj =" << NUM(totalConflicts) << "   best =" << NUM(bestSolutionValue)
	   << "   tenure =" << NUM(tabuTenure) 
	   << "   freq = " << NUMBER(frequency) 
	   << "   inc = " << NUM(increment)
	   << "   rand = " << NUM(randomTenure)
	   << "   max-min = " << NUM(max_min)
	   << endl;
      nextVerbose += incVerbose;
    }
  }

  freeArrays(nodesByColor, conflicts, tabuStatus, 
	     neighbors, nbcPosition, k, g.n);

  return result * totalIterations;
}




