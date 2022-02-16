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

#include "inputGraph.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <string>

using namespace std;

void inputDimacsGraph(Graph & g, char * file) {

    char c;
    char str[400];
    ifstream IN(file, ios::in);
    int line = 0;
    g.nbEdges = 0;
    int edges = -1;
    int blem = 1;
    int multiple = 0;

    std::string dump;

    while (!IN.eof()) {
        line++;
        IN.get(c);
        //    if (line < 100 || line % 1000 == 0) cout << "Read line " << line
        //    << " with code " << c <<  endl;
        if (IN.eof())
            break;

        switch (c) {
        case 'p':
            IN.get(c);
            IN.getline(str, 39, ' ');
            if (strcmp(str, "edge") && strcmp(str, "edges")) {
                cerr << "File " << file << " line " << line << ":\n";
                cerr << "Error reading 'p' line: no 'edge' keyword found.\n";
                cerr << "'" << str << "' found instead\n";
                exit(-1);
            }
            IN >> g.n;
            IN >> edges;
            blem = 0;
            g.resize(g.n);
            break;
        case 'n':
            if (blem) {
                cerr << "File " << file << " line " << line << ":\n";
                cerr << "Found 'n' line before a 'p' line.\n";
                exit(-1);
            }
            int node;
            IN >> node;
            IN >> dump;
            if (node < 1 || node > g.n) {
                cerr << "File " << file << " line " << line << ":\n";
                cerr << "Node number " << node << " is out of range!\n";
                exit(-1);
            }
            node--;
            // cout << "Tags (n Lines) not implemented in g object\n";
            // IN >> tag[node];
            break;
        case 'e':
            int node1, node2;
            IN >> node1 >> node2;
            if (node1 < 1 || node1 > g.n || node2 < 1 || node2 > g.n) {
                cerr << "File " << file << " line " << line << ":\n";
                cerr << "Node " << node1 << " or " << node2
                     << " is out of range!\n";
                exit(-1);
            }
            node1--;
            node2--;
            // Undirected graph
            // cout << "Read edge from " << node1 << " to " << node2 << endl;
            // if (g[node1][node2] == 0) {
            //     g.nbEdges++;
            // } else {
            //     multiple++;
            //     if (multiple < 5) {
            //         cerr << "Warning: in graph file " << file << " at line "
            //              << line << ": edge is defined more than once.\n";
            //         if (multiple == 4) {
            //             cerr << "  No more multiple edge warnings will be "
            //                     "issued\n";
            //         }
            //     }
            // }
						if(node1 != node2) {
							g.nbEdges++;

            	g.neighbor[node1].push_back(node2);
            	g.neighbor[node2].push_back(node1);
						}
            // g[node1][node2] = 1;
            // g[node2][node1] = 1;
            break;
        case 'd':
        case 'v':
        case 'x':
            cerr << "File " << file << " line " << line << ":\n";
            cerr << "'" << c << "' lines are not implemented yet...\n";
            IN.getline(str, 399, '\n');
            break;
        case 'c':
            IN.putback('c');
            IN.get(str, 399, '\n');
            break;
        default:
            cerr << "File " << file << " line " << line << ":\n";
            cerr << "'" << c << "' is an unknown line code\n";
            exit(-1);
        }
        IN.get(); // Kill the newline;
    }
    IN.close();
    if (multiple) {
        cerr << multiple << " multiple edges encountered\n";
    }

    int count = 0;
    for (int i = 0; i < g.n; ++i) {
        std::sort(g.neighbor[i].begin(), g.neighbor[i].end());
        count += g.neighbor[i].size();
        for (int j = 1; j < g.neighbor[i].size(); ++j) {
            assert(g.neighbor[i][j - 1] != g.neighbor[i][j]);
        }
    }

    std::cout << count / 2 << " / " << g.nbEdges << std::endl;
}

