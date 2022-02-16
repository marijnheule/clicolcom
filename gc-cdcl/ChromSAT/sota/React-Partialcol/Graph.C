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
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

Graph::Graph() {
    // matrix=NULL;
    n = 0;
    nbEdges = 0;
}

Graph::Graph(int m) {
    // matrix=NULL;
    resize(m);
}

// int * Graph::operator[](int index) {
//   if (index<0 || index >= this->n) {
//     cerr << "First node index out of range: " << index << "\n";
//     // Look for a proper way of treating this error...
//     matrix[-1]=0; //Make it crash.
//   }
//   return matrix+this->n*index;
// }

//  Graph::Graph(char * file) {
//    matrix=NULL;

//    this->n=0;
//    char c;
//    char str[400];
//    ifstream IN(file, ios::in);
//    int line=0;
//    this->nbEdges=0;
//    int edges=-1;
//    int blem=1;
//    while(!IN.eof()) {
//      line++;
//      IN.get(c);
//      if (IN.eof()) break;
//      switch (c) {
//      case 'p':
//        IN.get(c);
//        IN.getline(str,39,' ');
//        if (strcmp(str,"edge")) {
//  	cerr << "File " << file << " line " << line << ":\n";
//  	cerr << "Error reading 'p' line: no 'edge' keyword found.\n";
//  	cerr << "'" << str << "' found instead\n";
//  	exit(-1);
//        }
//        IN >> this->n;
//        IN >> edges;
//        cout << "Read " << this->n << " nodes and " << edges << " edges\n";
//        blem=0;
//        this->resize(this->n);
//        break;
//      case 'n':
//        if (blem) {
//  	cerr << "File " << file << " line " << line << ":\n";
//  	cerr << "Found 'n' line before a 'p' line.\n";
//  	exit(-1);
//        }
//        int node;
//        IN >> node;
//        if (node < 1 || node > this->n) {
//  	cerr << "File " << file << " line " << line << ":\n";
//  	cerr << "Node number " << node << " is out of range!\n";
//  	exit(-1);
//        }
//        node--;	
//        cout << "Tags (n Lines) not implemented in this object\n";
//        //IN >> tag[node];
//        break;
//      case 'e':
//        int node1, node2;
//        IN >> node1 >> node2;
//        if (node1 < 1 || node1 > this->n || node2 < 1 || node2 > this->n) {
//  	cerr << "File " << file << " line " << line << ":\n";
//  	cerr << "Node " << node1 << " or " << node2 << " is out of range!\n";
//  	exit(-1);
//        }
//        node1--;
//        node2--;
//        // Undirected graph
//        // cout << "Read edge from " << node1 << " to " << node2 << endl;
//        matrix[node1+this->n*node2] =1;
//        matrix[node2+this->n*node1] =1;
//        break;
//      case 'd':
//      case 'v':
//      case 'x':
//        cerr << "File " << file << " line " << line << ":\n";
//        cerr << "'" << c << "' lines are not implemented yet...\n";
//        IN.get(str,399,'\n');
//        break;
//      case 'c':
//        IN.get(str,399,'\n');
//        break;
//      default:
//  	cerr << "File " << file << " line " << line << ":\n";
//  	cerr << "'" << c << "' is an unknown line code\n";
//  	exit(-1);
//      }
//      IN.get(); // Kill the newline;
//    }
//    IN.close();
//  }

void Graph::resize(int m) {

    // std::cout << "resize " << m << std::endl;
    if (m > 0) {
        neighbor.clear();
        neighbor.resize(m);

        //   if (matrix != NULL) {
        //       delete[] matrix;
        // }

        n = m;
        nbEdges = 0;
        // matrix = new int[m*m];
        // for (int i=0; i<m*m; i++) {
        //   matrix[i]=0;
        // }
  }
}


Graph::~Graph() {
  resize(0);
}
