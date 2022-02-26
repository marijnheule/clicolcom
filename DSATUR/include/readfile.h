#ifndef __READFILE__
#define __READFILE__

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>

#include <vector>
#include <set>
#include <bitset>


using namespace std;

void usage();
void readFile(int &N, int &M, vector< vector<bool> > &edge );
void parseArgs(int argc, char **argv);

/***
 *  Procedures to relabel input graph (=DD variable ordering)
 ***/
bool RelabelGraphMinDist(vector< vector<int> >  &NeighborList,
			 vector< vector<bool> > &edge);
bool RelabelGraphDegree(vector< vector<int> >  &NeighborList,
		        vector< vector<bool> > &edge);
int RelabelDsatur(vector< vector<int> >  &NeighborList,
		  vector< vector<bool> > &edge,
		  bool onlyUpperBound, vector<int> &primalSol, vector<int> &vertexOrdering);


#endif
