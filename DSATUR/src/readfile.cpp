#include "../include/readfile.h"
#include <time.h>
#include <math.h>

extern const char* FileIn;
extern const char* FileOut; // to export solution
extern int timelimit;
extern int ordering;      // 0: lex order, 1: degree based (default), 2: distance based, 3: Dsatur
extern bool writesol;

#define DBG(x)

///////////////
void usage() {
  cout << endl
       << "USAGE: color -filein <filename>" << endl;
  exit(1);
}

///////////////////
void readFile(int &N, int &M, vector< vector<bool> > &edge ) {

  // dimaces format:
  //  c "random instance"
  //  p edges 120 3585
  //  e 1 7
  //  e 1 8
  //  e 1 10
  
  ifstream inFile(FileIn, ios::in);
  if (!inFile.is_open()) {
    cout << "error while opening file " << FileIn << endl;
    exit(1);
  }

  char chr;
  char str[255];

  N = -1; // Number of vertices
  M = -1; // This is the number of edges to be scanned in the file: may contain duplicates!

  int cnt = 0; // explicitly count unique edges
  
  inFile >> chr;
  while( chr=='c' || chr=='p' ) {
    if (chr=='c') {
      inFile.getline(str,2048);
      inFile >> chr;
    }
    else if (chr=='p') {
      inFile >> str; // read in "edges", "edge", or "col" (the input files are not consistent on this :(  )
      inFile >> N;
      //  cout << "Number of nodes: " << N << endl;
      inFile >> M;
      //  cout << "Number of edges: " << M << endl;
      inFile >> chr;
    }
  }

  if (N<=-1 || M<=-1) {
    cout << "error while reading file: no graph information in preamble" << endl;
  }
  edge.resize(N);
  for (int i=0; i<N; i++) {
    edge[i].resize(N);
  }

  for (int i=0; i<M; i++) {
    int e1, e2;
    if (i>=1) inFile >> chr; // "e"  (first 'e' was parsed already before)
    inFile >> e1;
    inFile >> e2;

    // check for duplicates
    if (edge[e1-1][e2-1] != 1) {    
      edge[e1-1][e2-1] = 1;
      edge[e2-1][e1-1] = 1;
      cnt++;
    }
  }    

  // return the true number of edges
  M = cnt;
  
  inFile.close();
}

////////////////////
void parseArgs(int argc, char **argv) {

  bool fileNameSpecified = false;

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    cout << argv[argIndex] << " ";
  }
  cout << endl;

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-filein") ) {
      argIndex++;
      if (argIndex > argc-1) {
	cout << "error: specify a value for filein" << endl;
	exit(1);
      }
      FileIn = argv[argIndex];
      fileNameSpecified = true;
    }
    else if ( !strcmp(argv[argIndex], "-solfile") ) {
      argIndex++;
      if (argIndex > argc-1) {
	cout << "error: specify a value for solfile" << endl;
	exit(1);
      }
      FileOut = argv[argIndex];
      writesol = true;
    }
    else {
      cout << "Unexpected option: " << argv[argIndex] << endl;
      usage();
    }
  }
  if (!fileNameSpecified) {
    cout << "please specify an input file" << endl;
    usage();
  }
}


/******************
 * relabel graph by reordering the nodes
 * based on minimum distance (bandwidth)
 ******************/


bool RelabelGraphMinDist(vector< vector<int> >  &NeighborList,
		  vector< vector<bool> > &edge) {
  /* ordering heuristic:
   * - evaluate pairwise distance for each edge
   * - apply 2-opt heuristic to minimize distance
   */

  //  cout << "relabeling graph based on minimum distance heuristic" << endl;

  bool success = false;
  int distance = 0;


//  vector<unsigned int> S;    // ordered set of selected nodes
//  set<unsigned int> Scomp;   // complement of S
//  vector<int> Deg;           // Deg[i] is degree of node i (not unsigned because need to compare to -1)
//  vector<int> Conn;          // Conn[i] is connectivity of i to S (number of neighbors in S)
  vector<unsigned int> Pos;  // Pos[i] is position of node i in the ordering

  // initialize
  int N = NeighborList.size();
  Pos.resize(N);
  for (int i=0; i<N; i++) {
    Pos[i] = i;
  }

  for (int i=0; i<N; i++) {
    for (std::vector<int>::iterator it = NeighborList[i].begin(); it != NeighborList[i].end(); ++it) {
      distance += abs((int)(Pos[i] - Pos[*it]));
    }
  }

  //cout << "  distance = " << distance << endl;
  
  bool cont = true;
  // apply 2-opt until no more improvement
  clock_t begin = clock();
  double elapsed_time = 0;
  
  while ((cont==true) && (elapsed_time <= 10)) {
    clock_t end = clock();
    elapsed_time = double(end - begin) / CLOCKS_PER_SEC;  


    int maxsav = 0; // max savings for pair of nodes (i,j)
    int maxi = -1;  // index of first node
    int maxj = -1;  // index of second node
    for (int i=0; i<N-1; i++) {
      for (int j=i+1; j<N; j++) {
        int sav_ij = 0; // local savings
        // impact of edges out of i
	for (std::vector<int>::iterator it = NeighborList[i].begin(); it != NeighborList[i].end(); ++it) {
          if ((*it != i) && (*it != j)) {
	    sav_ij += abs( (int)(Pos[*it] - Pos[i]) ) - abs( (int)(Pos[*it] - Pos[j]) );
	  }
	}
	// impact of edges out of j
	for (std::vector<int>::iterator it = NeighborList[j].begin(); it != NeighborList[j].end(); ++it) {
          if ((*it != i) && (*it != j)) {
            sav_ij += abs( (int)(Pos[*it] - Pos[j]) ) - abs( (int)(Pos[*it] - Pos[i]) );
	  }
	}
	if (sav_ij > maxsav) {
	  // update maximum 
	  maxsav = sav_ij;
	  maxi = i;
	  maxj = j;
	}
      }
    }
    if (maxsav <= 0) {
      //cout << "no further savings." << endl;
      cont = false;
    }
    else {
      // swap i and j
      //cout << "swap " << maxi << " and " << maxj << " with savings " << maxsav << endl;
      int tmp = Pos[maxi];
      Pos[maxi] = Pos[maxj];
      Pos[maxj] = tmp;
      
      //cout << "new ordering: ";
      //for (int i=0; i<N; i++) {
      //  cout << Pos[i] << " ";
      //}
      //cout << endl;

      // update distance to keep track
      distance -= maxsav;
      //cout << "distance updated to " << distance << endl;

      //cout << "  distance = " << distance << endl;

      if (distance < 0) {
        cout << "error in relabeling: negative total distance" << endl;
        exit (1);	
      }   
    }
  }

  // to relabel the graph, first copy the ordering
  // into the sequence S
  vector<unsigned int> S;    // ordered set of selected nodes
  S.resize(N);
  for (int i=0; i<N; i++) {
    S[Pos[i]] = i;
  }

  DBG(
  cout << "  finished relabeling with the following sequence:" << endl;
  for (std::vector<unsigned int>::iterator it = S.begin() ; it != S.end(); ++it)
    cout << " " << *it;
  cout << endl;
      )

  /* relabel the nodes and edges */
  vector< vector<int> > tmpNeighborList;
  tmpNeighborList.resize(N);
  for (int i=0; i<N; i++) {
    edge[i].clear();
    edge[i].resize(N);
    for (std::vector<int>::iterator it = NeighborList[S[i]].begin() ; it != NeighborList[S[i]].end(); ++it) {
      tmpNeighborList[i].push_back(Pos[*it]);
      edge[i][Pos[*it]] = 1;
    }
  }
  for (int i=0; i<N; i++) {
    NeighborList[i] = tmpNeighborList[i];
  }
 
  DBG( 
  cout << "edges after relabeling:" << endl;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      cout << " " << edge[i][j];
    }
    cout << endl;
  }
  cout << endl;
  )


  return success; 
}


/******************
 * relabel graph by reordering the nodes
 * based on degree and connectivity
 ******************/

bool RelabelGraphDegree(vector< vector<int> >  &NeighborList,
		        vector< vector<bool> > &edge) {
  /* ordering heuristic:
   *  - maintain ordered set S of nodes
   *  - among all neighbors of S select those with 
   *    * highest connectivity and
   *    * highest remaining degree
   */

  //  cout << "relabeling graph based on max-degree and connectivity" << endl;
  
  bool success = false;

  vector<unsigned int> S;    // ordered set of selected nodes
  set<unsigned int> Scomp;   // complement of S
  vector<int> Deg;           // Deg[i] is degree of node i (not unsigned because need to compare to -1)
  vector<int> Conn;          // Conn[i] is connectivity of i to S (number of neighbors in S)
  vector<unsigned int> Pos;  // Pos[i] is position of node i in sequence S

  // initialize
  int N = NeighborList.size();
  for (unsigned int i=0; i<N; i++) {
    Scomp.insert(i);
  }
  Deg.resize(N);
  for (int i=0; i<N; i++) {
    Deg[i] = NeighborList[i].size();
  }
  Conn.resize(N);
  Pos.resize(N);

  // add nodes to S using the static ordering heuristic described ahove
  for (int i=0; i<N; i++) {
    // select node with max connectivity and amx degree
    int maxConn = -1;         // maximum connectivity
    int maxDeg = -1;          // maximum degree
    int select = -1;          // index of selected node
    for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
      if (Conn[*it] > maxConn) {
        maxConn = Conn[*it];
        maxDeg = Deg[*it];
	select = *it;
      }
      else if ((Conn[*it] == maxConn) && (Deg[*it] > maxDeg)){
	maxDeg = Deg[*it];
        select = *it;
      }	
    }
    if (select == -1) {
      cout << "error in RelabelGraph for i = " << i << ": no node selected." << endl;
      exit(1);
    }

    // at this point, 'select' can be added to S and removed from Scomp
    S.push_back(select);
    Scomp.erase(select);
    // make updates to Deg, Conn, and Pos. 
    for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
      if ((edge[select][*it] == 1) || (edge[*it][select] == 1)) {
        Deg[*it] -= 1;
	Conn[*it] += 1;
      }
    }
    Pos[select] = S.size() - 1;
  }

  DBG(
  cout << "  finished relabeling with the following sequence:" << endl;
  for (std::vector<unsigned int>::iterator it = S.begin() ; it != S.end(); ++it)
    cout << " " << *it;
  cout << endl;
      )

  /* relabel the nodes and edges */
   

//bool RelabelGraph(vector< vector<int> >  &NeighborList,
//		  vector< vector<bool> > &edge) {
  vector< vector<int> > tmpNeighborList;
  tmpNeighborList.resize(N);
  for (int i=0; i<N; i++) {
    edge[i].clear();
    edge[i].resize(N);
    for (std::vector<int>::iterator it = NeighborList[S[i]].begin() ; it != NeighborList[S[i]].end(); ++it) {
      tmpNeighborList[i].push_back(Pos[*it]);
      edge[i][Pos[*it]] = 1;
    }
  }
  for (int i=0; i<N; i++) {
    NeighborList[i] = tmpNeighborList[i];
  }
 
  DBG( 
  cout << "edges after relabeling:" << endl;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      cout << " " << edge[i][j];
    }
    cout << endl;
  }
  cout << endl;
  )


  return success; 
}


/******************
 * relabel graph using Dsatur heuristic
 ******************/

int RelabelDsatur(vector< vector<int> >  &NeighborList,
		  vector< vector<bool> > &edge,
		  bool onlyUpperBound, vector<int> &primalSol, vector<int> &vertexOrdering) {
  /* ordering heuristic:
   *  - Use Dsatur heuristic [Brelaz, CACM 1979]
   *  - Select uncolored node with:
   *    o maximum saturated degree (number of distinct colors it is connected to)
   *    o in case of ties: select highest degree node
   *    o in case of ties: select highest connection to previously added nodes
   */

  //  cout << "relabeling graph based on Dsatur ordering" << endl;
  
  vector<unsigned int> S;    // ordered set of selected nodes
  set<unsigned int> Scomp;   // complement of S
  vector<int> Deg;           // Deg[i] is degree of node i (not unsigned because need to compare to -1)
  vector<unsigned int> Pos;  // Pos[i] is position of node i in sequence S
  vector<int> Conn;          // Conn[i] is connectivity of i to S (number of neighbors in S)

  vector<int> Col;           // Col[i] is the color assigned to node i (range is 1..N)

  vector< set<unsigned int> > Blocked; // Blocked[i] is the set of blocked colors for node i.
                                       // Use default ordering 'less' for sets (increasing).
  
  // initialize
  int N = NeighborList.size();
  for (unsigned int i=0; i<N; i++) {
    Scomp.insert(i);
  }
  Deg.resize(N);
  for (int i=0; i<N; i++) {
    Deg[i] = NeighborList[i].size();
  }
  Pos.resize(N);
  Conn.resize(N);
  Col.resize(N, 0);
  Blocked.resize(N);

  // Start by selecting a vertex of maximum degree
  int maxDeg = -1;          // maximum degree
  int select = -1;          // index of selected node
  for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
    if (Deg[*it] > maxDeg) {
      maxDeg = Deg[*it];
      select = *it;
    }
  }
  if (select == -1) {
    cout << "error in RelabelDsatur: no initial node selected." << endl;
    exit(1);
  }

  // add select to sequence, remove it from Scomp, and give it color 1
  S.push_back(select);
  Scomp.erase(select);
  Col[select] = 1;

  DBG(  cout << "select initial node " << select << " with degree " << maxDeg
	<< " and color it 1" << endl; )
  
  // make updates to Deg, Conn, Blocked, and Pos.
  for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
    if ((edge[select][*it] == 1) || (edge[*it][select] == 1)) {
      // 'select' is is a neighbor of '*it'.
      Deg[*it] -= 1;
      Conn[*it] += 1;
      Blocked[*it].insert(Col[select]);
    }
  }
  Pos[select] = S.size() - 1;
  
  // start adding the other vertices
  for (int i=1; i<N; i++) {
    // select node with max saturated degree, breaking ties with max connection and then max degree
    int maxSatDeg = -1;       // maximum saturated degree
    int maxDeg = -1;          // maximum degree
    int maxConn = -1;         // maximum connectivity
    int select = -1;          // index of selected node
    for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
      int size = Blocked[*it].size();
      if (size > maxSatDeg) {
        maxSatDeg = size;
        maxDeg = Deg[*it];
        maxConn = Conn[*it];
	select = *it;
      }
      else if ((size == maxSatDeg) && (Conn[*it] > maxConn)) {
        maxConn = Conn[*it];
	maxDeg = Deg[*it];
        select = *it;
      }	
      else if ((size == maxSatDeg) && (Conn[*it] == maxConn) && (Deg[*it] > maxDeg)) {
        maxConn = Conn[*it];
        select = *it;
      }	
    }
    if (select == -1) {
      cout << "error in RelabelGraph for i = " << i << ": no node selected." << endl;
      exit(1);
    }

    DBG(
    cout << "select next node " << select
	 << " with #Blocked = " << Blocked[select].size()
	 << ", Conn = " << Conn[select]
	 << ", Deg = " << Deg[select];
	)
    
    // at this point, 'select' can be added to S and removed from Scomp
    S.push_back(select);
    Scomp.erase(select);

    // give 'select' the lowest color that is not blocked
    for (int col=1; col<=N; col++) {
      if (Blocked[select].find(col) == Blocked[select].end()) {
	Col[select] = col;
	DBG(cout << " and color it " << Col[select] << endl;)
	break;
      }
    }
    
    // make updates to Deg, Conn, Blocked, and Pos.
    for (std::set<unsigned int>::iterator it=Scomp.begin(); it!=Scomp.end(); ++it) {
      if ((edge[select][*it] == 1) || (edge[*it][select] == 1)) {
	// 'select' is is a neighbor of '*it'.
	Deg[*it] -= 1;
	Conn[*it] += 1;
	Blocked[*it].insert(Col[select]);
      }
    }
    Pos[select] = S.size() - 1;
  }

  int MaxColor = 0;
  for (int i=0; i<N; i++) {
    if (Col[i] > MaxColor) {
      MaxColor = Col[i];
    }
  }
  cout << endl << "Dsatur upper bound: " << MaxColor << endl;

  // assign solution
  for (int i=0; i<N; i++) {
    primalSol[i] = Col[i];
  }
  
  
  /* relabel the nodes and edges */

  if (onlyUpperBound == 0) {

    DBG(
    cout << "  finished relabeling with the following sequence:" << endl;
    for (std::vector<unsigned int>::iterator it = S.begin() ; it != S.end(); ++it)
      cout << " " << *it;
    cout << endl;
	)

    // assign vertex ordering 
    int idx = 0;
    for (std::vector<unsigned int>::iterator it = S.begin() ; it != S.end(); ++it) {
      vertexOrdering[idx] = *it;
      idx++;
    }
      
    vector< vector<int> > tmpNeighborList;
    tmpNeighborList.resize(N);
    for (int i=0; i<N; i++) {
      edge[i].clear();
      edge[i].resize(N);
      for (std::vector<int>::iterator it = NeighborList[S[i]].begin() ; it != NeighborList[S[i]].end(); ++it) {
	tmpNeighborList[i].push_back(Pos[*it]);
	edge[i][Pos[*it]] = 1;
      }
    }
    for (int i=0; i<N; i++) {
      NeighborList[i] = tmpNeighborList[i];
    }
    
    DBG( 
	cout << "edges after relabeling:" << endl;
	for (int i=0; i<N; i++) {
	  for (int j=0; j<N; j++) {
	    cout << " " << edge[i][j];
	  }
	  cout << endl;
	}
	cout << endl;
	 )
  }

  return MaxColor;
}



