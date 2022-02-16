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

#include "Graph.h"
#include "Random.h"
#include "inputGraph.h"
#include "reactcol.h"
#include "verifyColoring.h"
#include "termcolors.h" 
#include "tabu.h"
#include <string.h>
#include <iostream>
#include <stdlib.h>


#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <limits.h>

using namespace std;

void usage() {
    cout << "Usage:\n";
    cout << "-k number of Colors\n"
         << "-f frequency\n"
         << "-i increment\n"
         << "-m maximal iterations without improvement\n"
         << "-r random seed\n"
         << "-T use the classical Tabucol algorithm instead of Reactcol\n"
         << "-t #tenure use a static tabu tenure instead of a reactive tenure\n"
         << "-x seconds (maximum allowed cpu time, default 8 hours)\n"
         << "[-v] (may be repeated to increase verbosity)\n"
         << "filname of graph to color";
}

int main(int argc, char ** argv) {

    Graph g;
    int k = 0;
    int frequency = 0;
    int increment = 0;
    int maxIterations = -1;
    int verbose = 0;
    int randomSeed = 43;
    int tenure = 0; // reactive tenure
    int algorithm = 1;
    int maxTime = 8 * 3600; // 8 hours cpu time max

    for (int i = 1; i < argc; i++) {

        if (strcmp("-k", argv[i]) == 0) {
            k = atoi(argv[++i]);

        } else if (strcmp("-f", argv[i]) == 0) {
            frequency = atoi(argv[++i]);

        } else if (strcmp("-i", argv[i]) == 0) {
            increment = atoi(argv[++i]);

        } else if (strcmp("-x", argv[i]) == 0) {
            maxTime = atoi(argv[++i]);

        } else if (strcmp("-m", argv[i]) == 0) {
            maxIterations = atoi(argv[++i]);

        } else if (strcmp("-t", argv[i]) == 0) {
            tenure = atoi(argv[++i]);

        } else if (strcmp("-T", argv[i]) == 0) {
            algorithm = 2;

        } else if (strcmp("-r", argv[i]) == 0) {
            randomSeed = atoi(argv[++i]);

        } else if (strcmp("-v", argv[i]) == 0) {
            verbose++;

        } else {
            inputDimacsGraph(g, argv[i]);
            cout << "Sucessfully read " << argv[i] << endl;
        }
    }

    if (k == 0) {
        usage();
        exit(-1);
    }

    if (maxIterations < 0) {
        maxIterations = g.n * g.n * 500;
    }

    int coloring[g.n];

    for (int i = 0; i < g.n; i++) {
        coloring[i] = 0;
    }

    int iter = 0;

    Random r(randomSeed);
    if (algorithm == 1) {
        iter = reactcol(g, coloring, k, r, maxIterations, tenure, verbose,
            frequency, increment, maxTime);
    } else {
        iter = tabu(g, coloring, k, r, maxIterations, tenure, verbose,
            frequency, increment, maxTime);
    }


		std::cout << "ITER=" << iter << std::endl;

    if (iter >= 0) {

        if (verbose) {
            if (!verifyColoring(g, coloring, k)) {
                cout << ERROR("This should not happen...") << endl;
                return -3;
            }
            cout << "Successful coloring with " << NUM(k) << " colors in "
                 << NUMBER(iter) << " iterations\n";
        } else {
            cout << "success after " << iter << " iterations\n";
        }
    } else {
        if (verbose)
            cout << "No coloring with " << NUM(k) << " colors found in "
                 << NUMBER(-iter) << " iterations\n";
        else
            cout << "Failed after " << -iter << " iterations\n";
    }
    if (!verbose) {
        cout << "colors " << k << endl;
    }

    // TIME IT
    tms myTimes; /* Variable qui va contenir les temps */
    long ticks;
    double sec;
    int min;
    times(&myTimes); /* Demande des valeurs */
    ticks = myTimes.tms_utime;
    sec = (double)ticks / (double)sysconf(_SC_CLK_TCK);
    min = (int)sec / 60;

    cout << "cpu_ms " << (int)(sec * 1000) << endl;
    cout << "maxIterations " << maxIterations << endl;
    cout << "tenure ";
    if (tenure)
        cout << tenure << endl;
    else
        cout << "reactive" << endl;

    return 0;
}



