#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <cstring>
#include <algorithm>
#include <math.h>

#include "../include/readfile.h"

#define DBG(x)
#define EPS 0.00001 // precision for LP-flow cut generation

/// global parameters
const char* FileIn;
const char* FileOut;
int timelimit = 3600;
int ordering = 1;        // 0: lex order, 1: degree based (default), 2: distance based, 3: Dsatur
int exact = 0;           // use exact (0: no (default), 1: exact, 2: restricted)

bool writesol = false;   // write best solution to file

using namespace std;

int main(int argc, char** argv) {

  DBG( cout << "Entering parseArgs..." << endl; )
  parseArgs(argc, argv);
  DBG( cout << "...leaving parseArgs" << endl; )

  cout << "Parameters:" << endl
       << "  filein = " << FileIn << endl
       << "  timelimit = " << timelimit << endl
       << "  ordering = " << ordering << endl;

  
  int N = 0; // number of graph nodes
  int M = 0; // number of graph edges
  vector< vector<bool> > edge;

  DBG( cout << "Entering readFile..." << endl; )
  readFile(N, M, edge);
  DBG( cout << "...leaving readFile" << endl; )

  cout << endl << "Input graph:" << endl
       << "  number of nodes: " << N << endl
       << "  number of edges: " << M << endl;

  vector< vector <int> > NeighborList;
  NeighborList.resize(N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      if (edge[i][j] == 1) {
	NeighborList[i].push_back(j);
      }
    }
  }
  
  clock_t begin = clock();
  double time_to_bestLB = -1;  // time to best found lower bound
  double time_to_bestUB = -1;  // time to best found upper bound

  int upper_bound = N;       // upper bound
  vector<int> primalSol(N);  // best primal solution (based on original graph labels), i.e., primalSol[i] is the color of vertex i
  for (int i=0; i<N; i++)
    primalSol[i] = 0;

  vector<int> vertexOrdering(N);  // vertex re-ordering
  for (int i=0; i<N; i++)
    vertexOrdering[i] = i; 

  upper_bound = RelabelDsatur(NeighborList, edge, 1, primalSol, vertexOrdering); // run Dsatur to get upper bound
  // maintain time to previous reporting event
  clock_t end = clock();
  double prev_time = double(end - begin) / CLOCKS_PER_SEC;  
  double elapsed_secs = prev_time;

  if (writesol) {
    ofstream solfile;
    solfile.open(FileOut, ios::out);
    for (int i=0; i<N; i++)
      solfile << primalSol[i] << endl;
    solfile.close();    
  }

  return 0;  
}




