/******************************************************************************/
/*                Graph coloring using techologies of Satz                    */
/*                Authors: Chumin LI & Quan Zhe                               */
/*                Copyright MIS, University of Picardie Jules Verne, France   */
/*                chu-min.li@u-picardie.fr                                    */
/*                Avril 2008                                                  */
/******************************************************************************/

/* The program is distributed for research purpose, but without any guarantee.
   For any other (commercial or industrial) use of the proram, please contact
   Chumin LI
*/

/* The program read a graph from a file in DIMACS format, a number of colors
   and the a heuristic
   Three heuristics for branching are available: dom, dom+degree, 
   dom+degree+UP
   Exploitation of symmetries between colors
*/

/* Based on color3.c, using parse parameters
 */

/* Based on color6, use node_state au lieu de nodeMark pour G_v
iMaxClique2 is working version in progress.
This is very different from iMaxClique1 which is a final version
 */

/* Based on iMaxClique2, simple lookahead after each backtracking
using matching (pair of nonneibor nodes)
*/

/* Based on iMaxClique3, lookahead after each backtracking
using isets using arrange_sets (insert nodes into an iset
in which a node is chosen as a candidate
Works well when the density is 90%
*/

/* Based on iMaxClique4, use iset_testable1 of iMaxClique1 to
control lookahead when choosing the next node to expand
Also use partitionIntoIsets() of iMaxClique1 when solving
a graph of low density
 */

/* Based on iMaxClique5, distinguish dense and sparse graphs
Works better for sparse graphs
 */

/* Based on iMaxClique7, compute node_nb_neibors after lookahead
also for dense graphs
 */

/* Based on iMaxClique8, when compute nb of neibors for dense graph
directly compute non neibors instead of from noneibors made passive
during expansion

The most efficient on the 22/09/09
 */

/* Based on iMaxClique14, use a new upper bound inspired from 
MaxSAT, limited to propagate unit isets in a partition for lookahead
*/

/* Based on ImaxClique18, propagate binary isets
 */

/* Based on ImaxClique19, apply more maxsatz
 */

/* Based on iMaxClique20, add a new node to every iset in an
inconsistent set.
 */

/* Based on iMaxClique21, test isets of size 2, 3, 4, 5 instead of 
only 2 in iMaxClique21
 */

/* Based on imaxClique22, when nb_extra_isets is one more than 
nb_conflicts, remove failed nodes and propagate unit iset
 */

/* Based on iMaxClique23, move lookahead from choose_candidate
to backtracking
*/

/* Based on iMaxClique23bis, when the graph is sparse, use
get_upper_bound_by_coloring
 */

/* Based on iMaxClique23bisbis, use adjacence matrice in 
patitionIntoIsetsTomita
*/

/* Based on iMaxClique31, expand clq according to adjacence matrice
instead of noneibors
*/

/* Based on iMaxClique32, compute degree information using adjacence
matrice for sparse graphs
 */

/* Based on iMaxClique33, use qsort in partitionIntoIsetTomita
 */

/* Based on iMaxClique35, simplify choose_candidate() function
 */

/* Based on iMaxClique36, clean the code cource
when partitioning the graph, insert the current node into
an iset of the maximum size
 */

/* Based on MaxCLQ1
 */

/* Based on MaxCLQ6, cleaning code
 */

/* Based on MaxCLQ8, nb_all_steps is initialized to NB_NODE instead of 0
to make the first branch computing node degree and sorting nodes 
 */

/* This is tryMaxCLQ8bis
 */

/* Based on MaxCLQ17, improve backtracking for maxclique
 */

/* Based on MaxCLQ18, introducing static_nb_noneibors[node], then in compute node
degree, compute node degree from candidate_stack or from non-neibors list, according
to the number of non neibors
*/

/* Based on MaxCLQ19, add extended failed literal detection
 */

/* Based on MaxCLQ23, when detect failed literals in an iset I, the detection should 
use in priority isets already used in previous failed literals in I (except 
when detecting the first failed literal in I
*/

/* Based on MaxCLQ24. (1) after partitioning into isets and running
MaxSatz, directly set candidate_iset_no when recopying vertices into
CLIQUE_CANDIDATE_STACK; (2) in filter(), search for isets of size 2 
in degree decreasing order of vertices (instead of increasing order 
as in MaxCLQ23 and MaxCLQ24
*/  

/* written from MaxCLQ25. Sort vertices in decreasing degree order. Search 
for cliques in subgraphs formed by vertices in this order by adding each time
one vertex (the subgraph is empty at the beginning). The exact maximum clique 
size or a upper bound is established for the subgraph containing the newly
added vertex. This upper bound is then used to prunning subtrees for subgraphes
containing this vertex and new added vertices. This idea is applied recursively,
i.e. in subgraphes containing branching vertcies
*/

/* Based on MaxCLQ25-cliquer, improved upper bound for a vertex, taking into
account not only the neighbors of the vertex, but also earlier UB
*/

/* Based on MaxCLQ25-cliquer2, use iset partition to see if the neibors of
the candidate can give a sufficiently large clique
*/

/* Based on MaxCLQ25-cliquer3, candidate nodes in the first 
MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer isets are directly get solved
 */

/* Based on MaxCLQ25-cliquer4, use maxsatz before expanding a node
 */

/* Based on MaxCLQ25-cliquer8, sort the input graph by increasing
degree order in successive subgraphs (insert the smallest degree node, 
then remove it from the input graph, insert the smallest degree
node, in the remaining graph, remove it, and so on...
*/

/* Based on MaxCLQ25-cliquer11 or MaxCLQ25-cliquer8, sort the input graph
by maximum iset (maximum clique in the complement input graph)
*/

/* Based on MaxCLQ25-cliquer11-MaxIsetOrder, sort each initial iset
in increasing degree order
*/ 

/* Based on MaxCLQ25-cliquer11-MaxIsetOrder2, incremental degree sorting 
of initial isets
*/ 

/* Based on MaxCLQ25-cliquer11-MaxIsetOrder3 and  MaxCLQ25-cliquer12
*/ 

/* Based on MaxCLQ25-cliquer12-MaxIsetOrderbis 
Further node test in  test_node_for_failed_nodes(int node, int iset)
*/ 

/* Based on MaxCLQ25-cliquer16-maxIsetOrder, in test_node_for_failed_nodes(int node, int iset)
use simple_further_test_node to avoid lookback and to simplify isets as soon as a failed node
is detected
*/

/* Based on MaxCLQ25-cliquer17-maxIsetOrder. For RB, do not seach for an initial maximal clique using a greedy heuristic
*/

/* Based on New2MaxCLQ25-cliquer18-maxIsetOrderSumDegree
 */

/* Based on IncMaxCLQ

- using qsort in build_complement_graph to sort neighbors and non-neighbors of a
  vertex
- when extra_iset_nb is 1, remove all vertices that have no support in an iset
  of size <=4
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>

#include <unistd.h>

#include <sys/resource.h>

#include <math.h>

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

#define WORD_LENGTH 100 
#define TRUE 1
#define FALSE 0
#define NONE -1
#define NONE2 -2

/* the tables of nodes and edges are statically allocated. Modify the 
   parameters tab_node_size and tab_edge_size before compilation if 
   necessary */
#define tab_node_size  10000
#define tab_edge_size 50000000
#define max_nb_value 100
#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item

#define PASSIVE 0
#define ACTIVE 1

int *node_neibors[tab_node_size];
int node_value[tab_node_size];
int node_state[2*tab_node_size];
//int *node_value_state[tab_node_size];
//int node_nb_value[tab_node_size];
//int node_nb_symmetric_value[tab_node_size];
//int node_impact[tab_node_size];
int nb_neibors[tab_node_size];
int static_nb_neibors[tab_node_size];
int static_nb_noneibors[tab_node_size];
//int branching_node[tab_node_size];
int NODE_STACK[tab_node_size];
int NODE_STACK_fill_pointer=0;

int UNIT_ISET_STACK[tab_node_size];
int UNIT_ISET_STACK_fill_pointer=0;

int edge[tab_edge_size][2];
int symmetry[max_nb_value];

int NB_VALUE, NB_NODE, NB_EDGE, H;

int NB_UNIT=0, NB_BACK = 0, NB_BRANCHE=0;

int non_neibor[tab_node_size];
int *node_non_neibors[tab_node_size];
int nb_non_neibor[tab_node_size];

int matrice[tab_node_size][tab_node_size];

int FORMAT=1;

my_type build_simple_graph_instance(char *input_file) {
  FILE* fp_in=fopen(input_file, "r");
  char ch, word2[WORD_LENGTH];
  int i, j, e, node1, node2, *neibors, neibor, redundant;
  if (fp_in == NULL) return FALSE;
  if (FORMAT==1) {
    fscanf(fp_in, "%c", &ch);
    while (ch!='p') {
      while (ch!='\n') fscanf(fp_in, "%c", &ch);  
      fscanf(fp_in, "%c", &ch);
    }
    fscanf(fp_in, "%s%d%d", word2, &NB_NODE, &NB_EDGE);
  }
  else    fscanf(fp_in, "%d%d",  &NB_NODE, &NB_EDGE);

  for(i=0; i<NB_EDGE; i++) {
    if (FORMAT==1)
      fscanf(fp_in, "%s%d%d", word2, &edge[i][0], &edge[i][1]);
    else fscanf(fp_in, "%d%d", &edge[i][0], &edge[i][1]);
    if (edge[i][0]==edge[i][1]) {
      printf("auto edge %d over %d\n", i--, NB_EDGE--);
    }
    else {
      if (edge[i][0]>edge[i][1]) {
	e=edge[i][1]; edge[i][1]=edge[i][0]; edge[i][0]=e;
      }
      nb_neibors[edge[i][0]]++;
      nb_neibors[edge[i][1]]++;
    }
  }
  fclose(fp_in);

  for(i=1; i<=NB_NODE; i++) 
    for(j=1; j<=NB_NODE; j++) {
      matrice[i][j]=FALSE;
      matrice[j][i]=FALSE;
    }

  for (i=1; i<=NB_NODE; i++) {
    static_nb_neibors[i]=nb_neibors[i];
    node_neibors[i]= (int *)malloc((nb_neibors[i]+1) * sizeof(int));
    node_neibors[i][nb_neibors[i]]=NONE;
    nb_neibors[i]=0;
    // node_nb_value[i]=NB_VALUE;
    // node_nb_symmetric_value[i]=NB_VALUE;
    node_state[i]=ACTIVE;
    // node_value_state[i]= (int *)malloc(NB_VALUE * sizeof(int));
    // for(j=0; j<NB_VALUE; j++)
    //  node_value_state[i][j]=ACTIVE;
  }
  for(i=0; i<NB_EDGE; i++) {
    node1=edge[i][0];    node2=edge[i][1]; 
    neibors=node_neibors[node1]; redundant=FALSE;
    for(j=0; j<nb_neibors[node1]; j++) 
      if (neibors[j]==node2) {
	printf("edge redundant %d \n", i);
	redundant=TRUE;
	break;
      }
    if (redundant==FALSE) {
      node_neibors[node1][nb_neibors[node1]++]=node2;
      node_neibors[node2][nb_neibors[node2]++]=node1;
      matrice[node1][node2]=TRUE;
      matrice[node2][node1]=TRUE;
    }
  }
  for (i=1; i<=NB_NODE; i++) {
    static_nb_neibors[i]=nb_neibors[i];
    static_nb_noneibors[i]=NB_NODE-nb_neibors[i]-1;
    node_neibors[i][nb_neibors[i]]=NONE;
  }
  return TRUE;
}

#define depth 1200

int CLIQUE_CANDIDATE_STACK_fill_pointer=0;
int CLIQUE_CANDIDATE_STACK[depth*tab_node_size];
int SOLVED_NODE_STACK_fill_pointer=0; 
int SOLVED_NODE_STACK[depth*tab_node_size];
int candidate_potentiel[depth*tab_node_size];
//int candidate_nb_neibors[depth*tab_node_size];
// int candidate_nb_noneibors[depth*tab_node_size];
// int candidate_position[depth*tab_node_size];
// int colorMark[max_nb_value];
int nodeMark[tab_node_size];
//int candidate_state[tab_node_size];
//int saved_candidate_color_nb[tab_node_size];
int node_nb_neibors[tab_node_size];
//int ISET_STACK_fill_pointer=0;
//int ISET_STACK[depth*tab_node_size];
//int ISETS_STACK_fill_pointer=0;
//int ISETS_STACK[tab_node_size];

int node_nb_noneibors[tab_node_size];
int NONEIBOR_STACK[tab_node_size];
int NONEIBOR_STACK_fill_pointer=0;
int CLIQUE_STACK[tab_node_size];
int CLIQUE_STACK_fill_pointer=0;
int MAXCLIQUE_STACK[tab_node_size];
int MAXCLIQUE_STACK_fill_pointer=0;
int saved_node_value[tab_node_size];
int node_iset_no[tab_node_size];
int candidate_iset_no[depth*tab_node_size];
int candidate_ub[depth*tab_node_size];
int UB=-1, LB=-1, NB_COLOR=-1;

int saved_candidate_stack[tab_node_size];
int saved_noneibor_stack[tab_node_size];
int degree_neibors[tab_node_size];
int MAX_CLQ_SIZE=0;
int NB_BACK_CLIQUE=0;

int verifyClique(int debut) {
  int i, node, *neibors, neibor;

  for(node=1; node<=NB_NODE; node++)
    nodeMark[node]=0;

  for(i=debut; i<CLIQUE_STACK_fill_pointer; i++) {
    node=CLIQUE_STACK[i];
    if (nodeMark[node]!=i-debut) {
      printf("erreur clique in node %d\n", node);
      return FALSE;
    }
    neibors=node_neibors[node];
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors))
      nodeMark[neibor]++;
  }
  return TRUE;
}

void sort_neibors_inc(int *neibors) {
  int *neibors1, *neibors2, neibor1, neibor2;
  neibors1=neibors;
  for(neibor1=*neibors1; neibor1!=NONE; neibor1=*(++neibors1)) {
    neibors2=neibors1+1;
    for(neibor2=*neibors2; neibor2!=NONE; neibor2=*(++neibors2)) {
      if (nb_neibors[neibor1]>nb_neibors[neibor2]) {
	*neibors1=neibor2;
	*neibors2=neibor1;
	neibor1=neibor2;
      }
    }
  }
}

void sort_neibors_dec(int *neibors) {
  int *neibors1, *neibors2, neibor1, neibor2;
  neibors1=neibors;
  for(neibor1=*neibors1; neibor1!=NONE; neibor1=*(++neibors1)) {
    neibors2=neibors1+1;
    for(neibor2=*neibors2; neibor2!=NONE; neibor2=*(++neibors2)) {
      if (nb_neibors[neibor1]<nb_neibors[neibor2]) {
	*neibors1=neibor2;
	*neibors2=neibor1;
	neibor1=neibor2;
      }
    }
  }
}

void sort_noneibors_dec(int *neibors) {
  int *neibors1, *neibors2, neibor1, neibor2;
  neibors1=neibors;
  for(neibor1=*neibors1; neibor1!=NONE; neibor1=*(++neibors1)) {
    neibors2=neibors1+1;
    for(neibor2=*neibors2; neibor2!=NONE; neibor2=*(++neibors2)) {
      if (nb_non_neibor[neibor1]<nb_non_neibor[neibor2]) {
	*neibors1=neibor2;
	*neibors2=neibor1;
	neibor1=neibor2;
      }
    }
  }
}

void sort_noneibors_inc(int *neibors) {
  int *neibors1, *neibors2, neibor1, neibor2;
  neibors1=neibors;
  for(neibor1=*neibors1; neibor1!=NONE; neibor1=*(++neibors1)) {
    neibors2=neibors1+1;
    for(neibor2=*neibors2; neibor2!=NONE; neibor2=*(++neibors2)) {
      if (nb_non_neibor[neibor1]>nb_non_neibor[neibor2]) {
	*neibors1=neibor2;
	*neibors2=neibor1;
	neibor1=neibor2;
      }
    }
  }
}

static int nondegreecmp_inc(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degneib1, degneib2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=nb_non_neibor[*node1]; degree2=nb_non_neibor[*node2];
  if (degree1<degree2)
    return -1;
  else if (degree1==degree2) {
      return 0;
  }
  else return 1;
}

static int static_degreecmp_dec(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degneib1, degneib2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=nb_neibors[*node1]; degree2=nb_neibors[*node2];
  if (degree1>degree2)
    return -1;
  else if (degree1==degree2) {
      return 0;
  }
  else return 1;
}

void build_complement_graph() {
  int node, node1, *neibors, neibor, nb;

  for (node=1; node<=NB_NODE; node++) 
    non_neibor[node]=TRUE;
  for (node=1; node<=NB_NODE; node++) {
    neibors=node_neibors[node]; nb=0;
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) {
      non_neibor[neibor]=FALSE;
      nb++;
    }
    non_neibor[node]=FALSE;
    node_non_neibors[node]= (int *)malloc((NB_NODE-nb) * sizeof(int));
    nb=0;
    for (node1=1; node1<=NB_NODE; node1++) {
      if (non_neibor[node1]==TRUE) {
	node_non_neibors[node][nb++]=node1;
      }
      else non_neibor[node1]=TRUE;
    }
    node_non_neibors[node][nb]=NONE;
    nb_non_neibor[node]=nb;
    if (nb+static_nb_neibors[node] !=NB_NODE-1)
      printf("erreur graphe by complementation\n");
  }
  for (node=1; node<=NB_NODE; node++) {
    qsort(node_non_neibors[node], nb_non_neibor[node], sizeof(int), 
	  nondegreecmp_inc);
    qsort(node_neibors[node], nb_neibors[node], sizeof(int), 
	  static_degreecmp_dec);
/*     sort_noneibors_inc(node_non_neibors[node]); */
/*     sort_neibors_dec(node_neibors[node]); */
  }
}

void complement_graph() {
  int node, nb, *neibors, i, j;

  for (node=1; node<=NB_NODE; node++) {
    neibors=node_neibors[node];
    node_neibors[node]=node_non_neibors[node];
    node_non_neibors[node]=neibors;
    nb=nb_neibors[node];
    if (nb != static_nb_neibors[node])
      printf("erreur graphe nb neibor\n");
    nb_neibors[node]=nb_non_neibor[node];
    static_nb_neibors[node]=nb_non_neibor[node];
    nb_non_neibor[node]=nb;
  }
  
/*   for (node=1; node<=NB_NODE; node++) { */
/*     sort_noneibors_inc(node_non_neibors[node]); */
/*     sort_neibors_dec(node_neibors[node]); */
/*   } */
  for(i=1; i<=NB_NODE; i++) 
    for(j=1; j<=NB_NODE; j++) {
      if (i!=j) {
	if (matrice[i][j]==TRUE)
	  matrice[i][j]=FALSE;
	else 
	  matrice[i][j]=TRUE;
      }
    }
}

void scanone(int argc, char *argv[], int i, int *varptr) {
  if (i>=argc || sscanf(argv[i],"%i",varptr)!=1){
    fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
    exit(-1);
  }
}

int HELP_FLAG=FALSE;
char *INPUT_FILE;
int LBflag=TRUE, UBflag=TRUE;

void parse_parameters(int argc,char *argv[]) {
  int i, temp, j;
  if (argc<2)
    HELP_FLAG=TRUE;
  else 
    for (i=1;i < argc;i++) {
      if (strcmp(argv[i],"-nbColors") == 0) 
	scanone(argc,argv,++i,&NB_VALUE);
      else if (strcmp(argv[i],"-f") == 0)
	scanone(argc,argv,++i,&FORMAT);
      else if (strcmp(argv[i],"-lb") == 0)
	scanone(argc,argv,++i,&LBflag);
     else if (strcmp(argv[i],"-ub") == 0)
	scanone(argc,argv,++i,&UBflag);
      else if (strcmp(argv[i],"-help") == 0)
	HELP_FLAG=TRUE;
      else   
	  INPUT_FILE=argv[i];
    }
}

char* filename(char* input) {
  char c, *input1;
  int  nb, nb1;
  input1=input; nb=0;
  for(c=*input1; c!='\0'; c=*(input1+=1)) 
    if (c=='/') nb++;
  input1=input; nb1=0;
  for(c=*input1; c!='\0'; c=*(input1+=1)) {
    if (c=='/') nb1++;
    if (nb==nb1)
      return input1+1;
  }
  return input1;
}

void compute_candidate_nb_neibors(int debut) {
  int i,  *neibors, neibor, node;
  for(i=debut; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i];
    neibors=node_neibors[node]; node_nb_neibors[node]=0;
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) 
      if (node_state[neibor]==TRUE) {
	node_nb_neibors[node]++;
      }
    // candidate_nb_neibors[i]=node_nb_neibors[node];
  }
}

void compute_degree_neibors(int debut) {
  int i, *neibors, neibor, node;
  for(i=debut; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i];
    neibors=node_neibors[node]; degree_neibors[node]=0;
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors))
      degree_neibors[node]=degree_neibors[node]+node_nb_neibors[neibor];
  }
}

int choose_min_candidate(int debut) {
  int i, chosen_node=NONE, node, min;
  min=NB_NODE;
  for(i=debut; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i];
    if (node_state[node]==ACTIVE && node_nb_neibors[node]<min) {
      min=node_nb_neibors[node];
      chosen_node=node;
    }
  }
  return chosen_node;
}

int choose_max_candidate(int debut) {
  int i, chosen_node=NONE, node, max;
  max=-1;
  for(i=debut; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i];
    if (node_state[node]==ACTIVE && node_nb_neibors[node]>max) {
      max=node_nb_neibors[node];
      chosen_node=node;
    }
  }
  return chosen_node;
}

#define debutGraph() CLIQUE_CANDIDATE_STACK_fill_pointer-\
candidate_potentiel[CLIQUE_CANDIDATE_STACK_fill_pointer-1]

int expand_clq_from_node1_old(int node) {
  int *neibors, neibor, saved_candidate_nb, i, nb, 
    *noneibors, noneibor, candidate, saved_noneibor_fp, debut;
  debut=debutGraph(); neibors=matrice[node];
  saved_candidate_nb=CLIQUE_CANDIDATE_STACK_fill_pointer; nb=0;
  for(i=debut; i<saved_candidate_nb; i++) {
    candidate=CLIQUE_CANDIDATE_STACK[i];
    if (neibors[candidate]==TRUE) {
      //node is no longer in the new graph
      //  candidate_nb_neibors[CLIQUE_CANDIDATE_STACK_fill_pointer]=
      //	candidate_nb_neibors[i];
      nb++;
      candidate_potentiel[CLIQUE_CANDIDATE_STACK_fill_pointer]=nb;
      push(candidate, CLIQUE_CANDIDATE_STACK);
    }
    else {
      node_state[candidate]=PASSIVE;
      push(candidate, NONEIBOR_STACK);
    }
  }
  return nb;
}

// Compute a maximal clique in the graph stored in 
// CLIQUE_CANDIDATE_STACK and put the vertices of the clique in 
// CLIQUE_STACK
void iCliqueCore(int debut) {
  int i, node, candidate, saved_candidate_stack_fill_pointer, debutGraph;

  // compute_candidate_nb_neibors(debut);
  saved_candidate_stack_fill_pointer=0;
  NONEIBOR_STACK_fill_pointer=0;
  while(saved_candidate_stack_fill_pointer<CLIQUE_CANDIDATE_STACK_fill_pointer) {
    candidate=choose_max_candidate(saved_candidate_stack_fill_pointer);
    if (candidate!=NONE) {
      debutGraph=saved_candidate_stack_fill_pointer;
      saved_candidate_stack_fill_pointer=CLIQUE_CANDIDATE_STACK_fill_pointer;
      expand_clq_from_node1_old(candidate);
      push(candidate, CLIQUE_STACK);
      node_state[candidate]=PASSIVE;
    }
  }
  for(i=0; i<NONEIBOR_STACK_fill_pointer; i++)
    node_state[NONEIBOR_STACK[i]]=ACTIVE;
  NONEIBOR_STACK_fill_pointer=0;
}

// Compute a clique and put the vertices of the clique in 
// CLIQUE_STACK from debut.
int iClique(int debut) {
  int i, node, candidate, saved_candidate_stack_fill_pointer, debutGraph;
  CLIQUE_STACK_fill_pointer=debut;
  CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for (node=1; node<=NB_NODE; node++) {
    if (node_state[node]==ACTIVE) {
      push(node, CLIQUE_CANDIDATE_STACK);
    }
  }
  compute_candidate_nb_neibors(debut);
  compute_degree_neibors(debut);
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    candidate_potentiel[i]=i+1;
    // candidate_nb_neibors[i]=node_nb_neibors[CLIQUE_CANDIDATE_STACK[i]];
  }
  iCliqueCore(debut);
  return CLIQUE_STACK_fill_pointer;
}

int CLIQUE_SIZES[tab_node_size];

int DENSITY, NB_ISETS_TEST=0;
double RATIO;
int ISET_NB=0;
int ISETS_SIZE[tab_node_size];
int ISETS[tab_node_size][tab_node_size];

static int degreecmp_static(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degreeNeibors1, degreeNeibors2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=node_nb_neibors[*node1]; degree2=node_nb_neibors[*node2];
  degreeNeibors1=degree_neibors[*node1];
  degreeNeibors2=degree_neibors[*node2];
  if (degree1>degree2)
    return -1;
  else if (degree1==degree2) {
    if (degreeNeibors1>degreeNeibors2)
      return -1;
    else if (degreeNeibors1==degreeNeibors2)
      return 0;
    else return 1;
  }
  else return 1;
}

static int degreecmp(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degreeNeibors1, degreeNeibors2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=node_nb_neibors[*node1]; degree2=node_nb_neibors[*node2];
  if (degree1>degree2)
    return -1;
  else if (degree1==degree2) {
      return 0;
  }
  else return 1;
}

int FILTER_STACK_fill_pointer=0;
int FILTER_STACK[tab_node_size];
int filter_state[tab_node_size];

#define NO_CONFLICT -2
int iset_state[tab_node_size];
int REDUCED_ISET_STACK[tab_node_size];
int REDUCED_ISET_STACK_fill_pointer=0;
int ENLARGED_ISET_STACK[tab_node_size];
int ENLARGED_ISET_STACK_fill_pointer=0;
int UNITISET_STACK[tab_node_size];
int UNITISET_STACK_fill_pointer=0;
int MY_UNITISET_STACK[tab_node_size];
int MY_UNITISET_STACK_fill_pointer=0;
int node_reason[tab_node_size];
#define NO_REASON -3
int iset_involved[tab_node_size];
int REASON_STACK[tab_node_size];
int REASON_STACK_fill_pointer=0;
int *CONFLICT_ISETS[3*tab_node_size];
int CONFLICT_ISET_STACK[tab_node_size];
int CONFLICT_ISET_STACK_fill_pointer=0;
int ADDED_NB_NODE;
//int conflict_iset[tab_node_size];
int PASSIVE_ISET_STACK[tab_node_size];
int PASSIVE_ISET_STACK_fill_pointer=0;

void assign_node_value(int node, int value, int reason) {
  node_value[node]=value;
  node_state[node]=PASSIVE;
  push(node, NODE_STACK);
  node_reason[node]=reason;
}

int exclude_noneibor(int noneibor, int reason_iset) {
  int my_iset;
  my_iset=node_iset_no[noneibor];
  assign_node_value(noneibor, FALSE, reason_iset);
  if (iset_state[my_iset]==ACTIVE) {
    ISETS_SIZE[my_iset]--;
    push(my_iset, REDUCED_ISET_STACK);
    if (ISETS_SIZE[my_iset]==1)
      push(my_iset, MY_UNITISET_STACK);
    else if (ISETS_SIZE[my_iset]==0)
      return my_iset;
  }
  return NO_CONFLICT;
}


int fix_node_for_iset(int node, int iset) {
  int *noneibors, noneibor, my_iset, i, *adjacence;
  iset_state[iset]=PASSIVE;
  push(iset, PASSIVE_ISET_STACK);
  assign_node_value(node, TRUE, iset);
  if (nb_non_neibor[node]>(FILTER_STACK_fill_pointer<<1)) {
    adjacence=matrice[node];
    for(i=0; i<FILTER_STACK_fill_pointer; i++) {
      noneibor=FILTER_STACK[i];
      if (node_state[noneibor]==ACTIVE && adjacence[noneibor] !=TRUE) {
	if ((my_iset=exclude_noneibor(noneibor, iset)) != NO_CONFLICT)
	  return my_iset;
      }
    }
  }
  else {
    noneibors=node_non_neibors[node];
    for(noneibor=*noneibors; noneibor!=NONE; noneibor=*(++noneibors))
      if (node_state[noneibor]==ACTIVE) {
	if ((my_iset=exclude_noneibor(noneibor, iset)) != NO_CONFLICT)
	  return my_iset;
      }
  }
  return NO_CONFLICT;
}

int fix_addedNode_for_iset(int node, int iset) {
  int *isets, c_iset;
  iset_state[iset]=PASSIVE;
  push(iset, PASSIVE_ISET_STACK);
  assign_node_value(node, FALSE, iset);
  isets=CONFLICT_ISETS[node]; 
  for(c_iset=*isets; c_iset!=NONE; c_iset=*(++isets))
    if (iset_state[c_iset]==ACTIVE) {
      ISETS_SIZE[c_iset]--;
      push(c_iset, REDUCED_ISET_STACK);
      if (ISETS_SIZE[c_iset]==1)
	push(c_iset, MY_UNITISET_STACK);
      else if (ISETS_SIZE[c_iset]==0)
	return c_iset;
    }
  return NO_CONFLICT;
}
/*
int fix_anyNode_for_iset(int node, int iset) {
  if (node>NB_NODE)
    return fix_addedNode_for_iset(node, iset);
  else return fix_node_for_iset(node, iset);
}
*/
int fix_unitIset(int iset) {
  int node, *nodes;

  nodes=ISETS[iset];
  for(node=*nodes; node!=NONE; node=*(++nodes)) 
    if (node_state[node]==ACTIVE) {
      if (node>NB_NODE)
	return fix_addedNode_for_iset(node, iset);
      else
	return fix_node_for_iset(node, iset);
    }
  printf("erreur unitIset\n");
  return NO_CONFLICT;
}

int unitIsetProcess() {
  int i, iset, j, my_iset;
  for(i=0; i<UNITISET_STACK_fill_pointer; i++) {
    iset=UNITISET_STACK[i];
    if (iset_state[iset]==ACTIVE && ISETS_SIZE[iset]==1) {
      MY_UNITISET_STACK_fill_pointer=0;
      if ((my_iset=fix_unitIset(iset))!=NO_CONFLICT)
	return my_iset;
      else 
	for(j=0; j<MY_UNITISET_STACK_fill_pointer; j++) {
	  iset=MY_UNITISET_STACK[j];
	  if (iset_state[iset]==ACTIVE)
	    if ((my_iset=fix_unitIset(iset))!=NO_CONFLICT)
	      return my_iset;
	}
    }
  }
  MY_UNITISET_STACK_fill_pointer=0;
  return NO_CONFLICT;
}

int ISET_STACK[tab_node_size], ISET_STACK_fill_pointer=0;

void lookback_for_maxsatz(int iset) {
  int i, *nodes, node, reason_iset;
  REASON_STACK_fill_pointer=0;
  push(iset, REASON_STACK); iset_involved[iset]=TRUE;
  for(i=0; i<REASON_STACK_fill_pointer; i++) {
    iset=REASON_STACK[i];
    nodes=ISETS[iset];
    for(node=*nodes; node!=NONE; node=*(++nodes))
      if (node_value[node]==FALSE && 
	  node_reason[node] !=NO_REASON && 
	  iset_involved[node_reason[node]]==FALSE) {
	reason_iset=node_reason[node];
	push(reason_iset, REASON_STACK);
	//	node_reason[node]=NO_REASON;
	iset_involved[reason_iset]=TRUE;
      }
  }
  for(i=0; i<REASON_STACK_fill_pointer; i++) 
    iset_involved[REASON_STACK[i]]=FALSE;
}

int INVOLVED_ISET_STACK[tab_node_size];
int INVOLVED_ISET_STACK_fill_pointer=0;
int iset_used[tab_node_size];

void enlarge_involved_iset() {
  int i, iset;
  
  CONFLICT_ISETS[ADDED_NB_NODE]=
    &CONFLICT_ISET_STACK[CONFLICT_ISET_STACK_fill_pointer];
  node_state[ADDED_NB_NODE]=ACTIVE;
  for(i=0; i<REASON_STACK_fill_pointer; i++) {
    iset=REASON_STACK[i];
    if (ISETS[iset][ISETS_SIZE[iset]] != NONE)
      printf("erreur conflict iset\n");
    ISETS[iset][ISETS_SIZE[iset]++]=ADDED_NB_NODE;
    ISETS[iset][ISETS_SIZE[iset]]=NONE;
    push(iset, CONFLICT_ISET_STACK);
    push(iset, ENLARGED_ISET_STACK);
  }
  push(NONE, CONFLICT_ISET_STACK);
  ADDED_NB_NODE++;
}

void enlarge_stored_involved_isets() {
  int i, iset;
  CONFLICT_ISETS[ADDED_NB_NODE]=
    &CONFLICT_ISET_STACK[CONFLICT_ISET_STACK_fill_pointer];
  node_state[ADDED_NB_NODE]=ACTIVE;
  for(i=0; i<INVOLVED_ISET_STACK_fill_pointer; i++) {
    iset=INVOLVED_ISET_STACK[i];
    if (iset_involved[iset]==FALSE) {
      iset_involved[iset]=TRUE;
      if (ISETS[iset][ISETS_SIZE[iset]] != NONE)
	printf("erreur conflict iset\n");
      ISETS[iset][ISETS_SIZE[iset]++]=ADDED_NB_NODE;
      ISETS[iset][ISETS_SIZE[iset]]=NONE;
      push(iset, CONFLICT_ISET_STACK);
      push(iset, ENLARGED_ISET_STACK);
    }
  }
  push(NONE, CONFLICT_ISET_STACK);
  ADDED_NB_NODE++;
  for(i=0; i<INVOLVED_ISET_STACK_fill_pointer; i++) {
    iset_involved[INVOLVED_ISET_STACK[i]]=FALSE;
    iset_used[INVOLVED_ISET_STACK[i]]=FALSE;
  }
  INVOLVED_ISET_STACK_fill_pointer=0;
}

void reset_context_for_maxsatz(int saved_node_stack_fill_pointer,
			       int saved_passive_iset_stack_fill_pointer,
			       int saved_reduced_iset_stack_fill_pointer,
			       int saved_unitiset_stack_fill_pointer) {
  int i, node;

  for(i=saved_node_stack_fill_pointer; i<NODE_STACK_fill_pointer; i++) {
    node=NODE_STACK[i];
    node_state[node]=ACTIVE;
    node_reason[node]=NO_REASON;
  }
  NODE_STACK_fill_pointer=saved_node_stack_fill_pointer;
  for(i=saved_passive_iset_stack_fill_pointer; 
      i<PASSIVE_ISET_STACK_fill_pointer; i++)
    iset_state[PASSIVE_ISET_STACK[i]]=ACTIVE;
  PASSIVE_ISET_STACK_fill_pointer=saved_passive_iset_stack_fill_pointer;
  for(i=saved_reduced_iset_stack_fill_pointer; 
      i<REDUCED_ISET_STACK_fill_pointer; i++)
    ISETS_SIZE[REDUCED_ISET_STACK[i]]++;
  REDUCED_ISET_STACK_fill_pointer=saved_reduced_iset_stack_fill_pointer;
  UNITISET_STACK_fill_pointer=saved_unitiset_stack_fill_pointer;
  MY_UNITISET_STACK_fill_pointer=0;
}

int TESTED_NODE_STACK[tab_node_size];
int TESTED_NODE_STACK_fill_pointer=0;
int node_tested_state[tab_node_size];

void reset_context_for_maxsatz_no(int saved_node_stack_fill_pointer,
				  int saved_passive_iset_stack_fill_pointer,
				  int saved_reduced_iset_stack_fill_pointer,
				  int saved_unitiset_stack_fill_pointer) {
  int i, node;

  for(i=saved_node_stack_fill_pointer; i<NODE_STACK_fill_pointer; i++) {
    node=NODE_STACK[i];
    node_state[node]=ACTIVE;
    node_reason[node]=NO_REASON;
    if (node_value[node]==TRUE && node_tested_state[node]==FALSE) {
      node_tested_state[node]=TRUE; // no need to re-test at this point
      push(node, TESTED_NODE_STACK);
    }
  }
  NODE_STACK_fill_pointer=saved_node_stack_fill_pointer;
  for(i=saved_passive_iset_stack_fill_pointer; 
      i<PASSIVE_ISET_STACK_fill_pointer; i++)
    iset_state[PASSIVE_ISET_STACK[i]]=ACTIVE;
  PASSIVE_ISET_STACK_fill_pointer=saved_passive_iset_stack_fill_pointer;
  for(i=saved_reduced_iset_stack_fill_pointer; 
      i<REDUCED_ISET_STACK_fill_pointer; i++)
    ISETS_SIZE[REDUCED_ISET_STACK[i]]++;
  REDUCED_ISET_STACK_fill_pointer=saved_reduced_iset_stack_fill_pointer;
  UNITISET_STACK_fill_pointer=saved_unitiset_stack_fill_pointer;
  MY_UNITISET_STACK_fill_pointer=0;
}

void reset_enlarged_isets() {
  int i, iset;
  for(i=0; i<ENLARGED_ISET_STACK_fill_pointer; i++) {
    iset=ENLARGED_ISET_STACK[i];
    iset_state[iset]=ACTIVE;
    ISETS_SIZE[iset]--;
    ISETS[iset][ISETS_SIZE[iset]]=NONE;
  }
  ENLARGED_ISET_STACK_fill_pointer=0;
}

int my_unitIsetProcess() {
  int j, iset, iset_start=0, used_iset_start=0, my_iset;
  do {
    for(j=used_iset_start; j<MY_UNITISET_STACK_fill_pointer; j++) {
      iset=MY_UNITISET_STACK[j];
      if (iset_state[iset] == ACTIVE && iset_used[iset]==TRUE)
	if ((my_iset=fix_unitIset(iset))!=NO_CONFLICT)
	  return my_iset;
    }
    used_iset_start=j;
    for(j=iset_start; j<MY_UNITISET_STACK_fill_pointer; j++) {
      iset=MY_UNITISET_STACK[j];
      if (iset_state[iset] == ACTIVE) {
	if ((my_iset=fix_unitIset(iset))!=NO_CONFLICT)
	  return my_iset;
	iset_start=j+1;
	break;
      }
    }
  }
  while (j<MY_UNITISET_STACK_fill_pointer);
  return NO_CONFLICT;
}

int fix_anyNode_for_iset(int node, int iset) {
  if (node>NB_NODE)
    return fix_addedNode_for_iset(node, iset);
  else
    return fix_node_for_iset(node, iset);
}

void store_involved_isets() {
  int i, iset;
  for(i=0; i<REASON_STACK_fill_pointer; i++) {
    iset=REASON_STACK[i];
    push(iset, INVOLVED_ISET_STACK);
    iset_used[iset]=TRUE;
  }
}

int choose_iset(int start) {
  int i, chosen_iset_index=NONE, iset, chosen_iset=NONE, 
    min=2000000000, active_nodes[2], node, *nodes, ctr;
  for(i=start; i<REDUCED_ISET_STACK_fill_pointer; i++) {
    iset=REDUCED_ISET_STACK[i];
    if (iset_state[iset]==ACTIVE && ISETS_SIZE[iset]==2) {
      nodes=ISETS[iset]; ctr=0;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node_state[node]==ACTIVE) {
	  if (node<=NB_NODE) {
	    node_reason[node]=NO_REASON;
	    active_nodes[ctr++]=node;
	  }
	  else break;
	}
      }
      if (ctr==2 && node==NONE) 
	if (node_nb_neibors[active_nodes[0]]*node_nb_neibors[active_nodes[1]]<min) {
	  chosen_iset=iset; 
	  min=node_nb_neibors[active_nodes[0]]*node_nb_neibors[active_nodes[1]];
	}
    }
  }
  return chosen_iset;
}

int further_test_node(int start) {
  int my_iset, saved_unitiset_stack_fill_pointer, 
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer, chosen_iset, node, *nodes, i,
    iset_tested[tab_node_size];

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  for(i=start; i<saved_reduced_iset_stack_fill_pointer; i++)
    iset_tested[REDUCED_ISET_STACK[i]]=FALSE;
  for(i=start; i<saved_reduced_iset_stack_fill_pointer; i++) {
    chosen_iset=REDUCED_ISET_STACK[i];
    if (iset_tested[chosen_iset]==FALSE &&
	iset_state[chosen_iset]==ACTIVE && ISETS_SIZE[chosen_iset]==2) {
      nodes=ISETS[chosen_iset]; iset_tested[chosen_iset]=TRUE;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node_state[node]==ACTIVE && node>NB_NODE) 
	  break;
      }
      if (node==NONE) { 
	//chosen_iset=choose_iset(start);
	// if (chosen_iset != NONE) {
	nodes=ISETS[chosen_iset];
	for(node=*nodes; node!=NONE; node=*(++nodes)) {
	  if (node_state[node]==ACTIVE) {
	    MY_UNITISET_STACK_fill_pointer=0;
	    if ((my_iset=fix_node_for_iset(node, chosen_iset)) != NO_CONFLICT ||
		(my_iset=my_unitIsetProcess()) !=  NO_CONFLICT) {
	      iset_involved[chosen_iset]=TRUE;
	      lookback_for_maxsatz(my_iset);
	      iset_involved[chosen_iset]=FALSE;
	      reset_context_for_maxsatz(saved_node_stack_fill_pointer,
					saved_passive_iset_stack_fill_pointer,
					saved_reduced_iset_stack_fill_pointer,
					saved_unitiset_stack_fill_pointer);
	      store_involved_isets();
	    }
	    else {
	      reset_context_for_maxsatz(saved_node_stack_fill_pointer,
					saved_passive_iset_stack_fill_pointer,
					saved_reduced_iset_stack_fill_pointer,
					saved_unitiset_stack_fill_pointer);
	      break;
	    }
	  }
	}
	if (node==NONE)
	  return chosen_iset;
      }
    }
  }
  return  NO_CONFLICT;
}

int test_node(int node, int iset, int further) {
  int my_iset, saved_unitiset_stack_fill_pointer, 
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer;

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  MY_UNITISET_STACK_fill_pointer=0;
  //  if (ISETS_SIZE[iset]>1 && (my_iset=unitIsetProcess())!=NO_CONFLICT)
  //  printf("bizzar.....\n");
  if ((my_iset=fix_node_for_iset(node, iset)) != NO_CONFLICT ||
      (my_iset=my_unitIsetProcess()) !=  NO_CONFLICT ||
      (further &&
       (my_iset=further_test_node(saved_reduced_iset_stack_fill_pointer)) 
       != NO_CONFLICT)) {
    lookback_for_maxsatz(my_iset);
    reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			      saved_passive_iset_stack_fill_pointer,
			      saved_reduced_iset_stack_fill_pointer,
			      saved_unitiset_stack_fill_pointer);
    return my_iset;
  }
  else {
    reset_context_for_maxsatz_no(saved_node_stack_fill_pointer,
				 saved_passive_iset_stack_fill_pointer,
				 saved_reduced_iset_stack_fill_pointer,
				 saved_unitiset_stack_fill_pointer);
    return NO_CONFLICT;
  }
}

#define FL_TEST_LENGTH 3

int maxsatz_lookahead_by_fl(int nb_conflict, int nb_extra_isets) {
  int iset, *nodes, node, no_conflict, test_flag, k, saved_size, ctr, i;
  for(i=0; i<TESTED_NODE_STACK_fill_pointer; i++) 
    node_tested_state[TESTED_NODE_STACK[i]]=FALSE;
  TESTED_NODE_STACK_fill_pointer=0; 
  for(k=2; k<=FL_TEST_LENGTH; k++) 
    //for(i=0; i<=ISET_NB-nb_extra_isets+nb_conflict; i++) {
   for(iset=ISET_NB-1; iset>=nb_extra_isets-nb_conflict; iset--) {
    if (iset_state[iset]==ACTIVE && ISETS_SIZE[iset]==k) {
      nodes=ISETS[iset]; test_flag=TRUE;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node>NB_NODE || node_tested_state[node]==TRUE) {
	  test_flag=FALSE;
	  break;
	}
      }
      if (test_flag==TRUE) {
	nodes=ISETS[iset]; no_conflict=FALSE; 
	ctr=k+1;
	INVOLVED_ISET_STACK_fill_pointer=0;
/* 	if (INVOLVED_ISET_STACK_fill_pointer !=0) */
/* 	  printf("erreur involved_iset_stack"); */
	for(node=*nodes; node!=NONE; node=*(++nodes)) {
	  if (node_tested_state[node]==TRUE ||
	      test_node(node, iset, (*(nodes+1)==NONE))==NO_CONFLICT) {
	    no_conflict=TRUE;
	    break;
	  }
	  else store_involved_isets();
	}
	if (no_conflict==FALSE) { 
	  saved_size=ISETS_SIZE[iset];
	  enlarge_stored_involved_isets();
	  if (saved_size+1 != ISETS_SIZE[iset])
	    printf("erreur iset involved %d %d %d...", 
		   saved_size,  ISETS_SIZE[iset], iset);
	  nb_conflict++;
	  if (nb_extra_isets<=nb_conflict) {
/* 	    for(i=0; i<TESTED_NODE_STACK_fill_pointer; i++) */
/* 	      node_tested_state[TESTED_NODE_STACK[i]]=FALSE; */
/* 	    TESTED_NODE_STACK_fill_pointer=0; */
	    return nb_conflict;
	  }
	}
	for(i=0; i<INVOLVED_ISET_STACK_fill_pointer; i++)
	  iset_used[INVOLVED_ISET_STACK[i]]=FALSE;
	INVOLVED_ISET_STACK_fill_pointer=0;
      }
    }
  }
/*   for(i=0; i<TESTED_NODE_STACK_fill_pointer; i++) */
/*     node_tested_state[TESTED_NODE_STACK[i]]=FALSE; */
/*   TESTED_NODE_STACK_fill_pointer=0; */
  return nb_conflict;
}

static int simplecmp_dec(const void *p_int1, const void *p_int2) {
  int *pint1, *pint2, int1, int2;
  pint1=(int *) p_int1; pint2=(int *) p_int2;
  int1=*pint1; int2=*pint2;
  if (int1>int2)
    return 1;
  else if (int1==int2) {
      return 0;
  }
  else return -1;
}

int simple_further_test_node(int start) {
  int my_iset, saved_unitiset_stack_fill_pointer, 
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer, chosen_iset, node, *nodes, i, j;
  int my_saved_unitiset_stack_fill_pointer,  iset_tested[tab_node_size],
    my_saved_node_stack_fill_pointer, my_saved_passive_iset_stack_fill_pointer,
    my_saved_reduced_iset_stack_fill_pointer, conflict=FALSE, test_flag;

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  my_saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  my_saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  my_saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  my_saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
/*qsort(&(REDUCED_ISET_STACK[start]), REDUCED_ISET_STACK_fill_pointer-start, */
/* 	sizeof(int), simplecmp_dec); */
  for(i=start; i<REDUCED_ISET_STACK_fill_pointer; i++)
    iset_tested[REDUCED_ISET_STACK[i]]=FALSE;
  //  do { test_flag=0;
  for(i=start; i<REDUCED_ISET_STACK_fill_pointer; i++) {
    chosen_iset=REDUCED_ISET_STACK[i];
    if (iset_state[chosen_iset]==ACTIVE && iset_tested[chosen_iset]==FALSE &&
	ISETS_SIZE[chosen_iset]<=2) {
      nodes=ISETS[chosen_iset]; iset_tested[chosen_iset]=TRUE;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node<=NB_NODE && node_state[node]==ACTIVE) {
	  MY_UNITISET_STACK_fill_pointer=0;
	  my_iset=fix_node_for_iset(node, chosen_iset);
	  if (my_iset == NO_CONFLICT)
	    my_iset=my_unitIsetProcess();
	  reset_context_for_maxsatz(my_saved_node_stack_fill_pointer,
				    my_saved_passive_iset_stack_fill_pointer,
				    my_saved_reduced_iset_stack_fill_pointer,
				    my_saved_unitiset_stack_fill_pointer);
	  if (my_iset != NO_CONFLICT) {
	    assign_node_value(node, FALSE, NO_REASON);
	    ISETS_SIZE[chosen_iset]--;
	    push(chosen_iset, REDUCED_ISET_STACK);
	    if (ISETS_SIZE[chosen_iset]==1) {
	      MY_UNITISET_STACK_fill_pointer=0;  test_flag++;
	      push(chosen_iset, MY_UNITISET_STACK);
	      if (my_unitIsetProcess() != NO_CONFLICT) {
		conflict=TRUE;
		break;
	      }
	      for(j=my_saved_reduced_iset_stack_fill_pointer;
		  j<REDUCED_ISET_STACK_fill_pointer; j++)
		iset_tested[REDUCED_ISET_STACK[j]]=FALSE;
	      my_saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
	      my_saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
	      my_saved_passive_iset_stack_fill_pointer=
		PASSIVE_ISET_STACK_fill_pointer;
	      my_saved_reduced_iset_stack_fill_pointer=
		REDUCED_ISET_STACK_fill_pointer;
	    }
	  }
	}
      }
      if (conflict==TRUE)
	break;
    }
  }
  //  } while (conflict==FALSE && test_flag>0);
  reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			    saved_passive_iset_stack_fill_pointer,
			    saved_reduced_iset_stack_fill_pointer,
			    saved_unitiset_stack_fill_pointer);
  if (conflict==TRUE)
    return chosen_iset;
  else return NO_CONFLICT;
}

int test_node_for_failed_nodes(int node, int iset) {
  int my_iset, saved_unitiset_stack_fill_pointer, 
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer;

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  MY_UNITISET_STACK_fill_pointer=0;
  if ((my_iset=fix_node_for_iset(node, iset)) == NO_CONFLICT) {
    if ((my_iset=my_unitIsetProcess()) == NO_CONFLICT)
      my_iset=simple_further_test_node(saved_reduced_iset_stack_fill_pointer);
  }
  reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			    saved_passive_iset_stack_fill_pointer,
			    saved_reduced_iset_stack_fill_pointer,
			    saved_unitiset_stack_fill_pointer);
  return my_iset;
}

void check_consistency() {
  int iset, node, *nodes, nb_a, nb_p, nb, nb_passive_nodes=0;

  nb_p=0;
  for(iset=0; iset<ISET_NB; iset++) {
    if (iset_involved[iset] !=FALSE)
      printf("erreur involved... ");
    if (iset_state[iset]==ACTIVE) {
      nodes=ISETS[iset]; nb_a=0; 
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node_state[node]==ACTIVE) 
	  nb_a++;
	else if (node<=NB_NODE)
	  nb_passive_nodes++;
      }
      if (ISETS_SIZE[iset]!=nb_a)
	printf("erreur nb_a");
    }
    else {
      nb_p++;
      nodes=ISETS[iset]; nb_a=0; nb=0;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node_state[node]!=ACTIVE && node<=NB_NODE)
	  nb_passive_nodes++;
	if (node_state[node]==ACTIVE) 
	  printf("erreur active...");
	else if (node<=NB_NODE && node_value[node]==TRUE)
	  nb_a++;
	else if (node>NB_NODE)
	  nb++;
      }
      if (nb_a==0 && nb==0)
	printf("erreur SAT...");
    }
  }
  // if (nb_p != REDUCED_ISET_STACK_fill_pointer)
  //  printf("erreur nb_p %d %d\n", nb_p, REDUCED_ISET_STACK_fill_pointer);
  for(node=NB_NODE+1; node<ADDED_NB_NODE; node++)
    if (node_state[node]==PASSIVE)
      nb_passive_nodes++;
  if (nb_passive_nodes != NODE_STACK_fill_pointer)
    printf("erreur active node...");
  if (nb_p != PASSIVE_ISET_STACK_fill_pointer)
    printf("erreur nb_p %d %d\n", nb_p, PASSIVE_ISET_STACK_fill_pointer);
}

int test_by_eliminate_failed_nodes() {
  int node, my_iset, *nodes, saved_unitiset_stack_fill_pointer,
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer, conflict, false_flag;
  // return TRUE;
  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  do {false_flag=0;
  for(my_iset=ISET_NB-1; my_iset>=0; my_iset--) {
    // for(my_iset=0; my_iset<ISET_NB; my_iset++) {
    if (iset_state[my_iset]==ACTIVE) {
      //  check_consistency();
      nodes=ISETS[my_iset]; conflict=FALSE; 	 
      MY_UNITISET_STACK_fill_pointer=0;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node<=NB_NODE && node_state[node]==ACTIVE && 
	    test_node_for_failed_nodes(node, my_iset)!=NO_CONFLICT) {
	  MY_UNITISET_STACK_fill_pointer=0;
	  assign_node_value(node, FALSE, NO_REASON);
	  ISETS_SIZE[my_iset]--; false_flag++;
	  push(my_iset, REDUCED_ISET_STACK);
	  if (ISETS_SIZE[my_iset]==1) {
	    push(my_iset, MY_UNITISET_STACK);
	    break;
	  }
	  else if (ISETS_SIZE[my_iset]==0) {
	    conflict=TRUE;
	    break;
	  }
	}
      }
      if (conflict==TRUE)
	break;
      else if (MY_UNITISET_STACK_fill_pointer>0 &&
	       my_unitIsetProcess() != NO_CONFLICT) {
	    conflict=TRUE;
	    break;
      }
      //      if (false_flag==TRUE) my_iset=0;
    }
  }
  } while (false_flag>1 && conflict==FALSE);
  reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			    saved_passive_iset_stack_fill_pointer,
			    saved_reduced_iset_stack_fill_pointer,
			    saved_unitiset_stack_fill_pointer);
  if (conflict==TRUE)
    return NONE;
  else
    return TRUE;
}

int NB1=0;

int maxsatz(int clq_size) {
  int i, nb_conflict=0, iset, nb=0, nb_node=0,
    saved_unitiset_stack_fill_pointer, nb_extra_isets, max_size,
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer;
  // return TRUE;
  if (clq_size <=0)
    return 0;
  UNITISET_STACK_fill_pointer=0; nb_extra_isets=ISET_NB-clq_size;
  max_size=0;
  for(i=0; i<ISET_NB; i++) {
    iset_used[i]=FALSE;
    if (max_size<ISETS_SIZE[i])
      max_size=ISETS_SIZE[i];
    iset_state[i]=ACTIVE;
    if (ISETS_SIZE[i]==1)
      push(i, UNITISET_STACK);
    else if  (ISETS_SIZE[i]==2)
      nb++;
  }
  // if (nb_node != candidate_potentiel[CLIQUE_CANDIDATE_STACK_fill_pointer-1])
  // printf("erreur isets_size\n");
  // printf("%d %d %d %d\n", 
  //	  ISET_NB, clq_size, UNITISET_STACK_fill_pointer, nb);
  
  // if (UNITISET_STACK_fill_pointer/2+nb/5<nb_extra_isets)
  //  return TRUE;
  ADDED_NB_NODE=NB_NODE+1;
  CONFLICT_ISET_STACK_fill_pointer=0;
  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  NODE_STACK_fill_pointer=0; saved_node_stack_fill_pointer=0;
  PASSIVE_ISET_STACK_fill_pointer=0;
  saved_passive_iset_stack_fill_pointer=0;
  REDUCED_ISET_STACK_fill_pointer=0;
  saved_reduced_iset_stack_fill_pointer=0;
  ENLARGED_ISET_STACK_fill_pointer=0;
  MY_UNITISET_STACK_fill_pointer=0;
  
  while ((iset=unitIsetProcess())!=NO_CONFLICT) {
    lookback_for_maxsatz(iset);
    reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			      saved_passive_iset_stack_fill_pointer,
			      saved_reduced_iset_stack_fill_pointer,
			      saved_unitiset_stack_fill_pointer);
    enlarge_involved_iset();
    nb_conflict++;
    if (nb_extra_isets<=nb_conflict)
      break;
  }
  reset_context_for_maxsatz_no(saved_node_stack_fill_pointer,
			       saved_passive_iset_stack_fill_pointer,
			       saved_reduced_iset_stack_fill_pointer,
			       saved_unitiset_stack_fill_pointer);
  //  if (nb_extra_isets>nb_conflict)   
  //  nb_conflict=maxsatz_lookahead_by_fl_for_unitiset(nb_conflict, nb_extra_isets);
  if (nb_extra_isets>nb_conflict) 
    nb_conflict=maxsatz_lookahead_by_fl(nb_conflict, nb_extra_isets);
  // for(i=0; i<ISET_NB; i++) 
  //  nb_node -= ISETS_SIZE[i];
  // if (nb_node != 0)
  //  printf("erreur isets_size -\n");
  if (nb_conflict != ADDED_NB_NODE-NB_NODE-1)
    printf("erreur nb conflict %d %d...\n", nb_conflict,ADDED_NB_NODE-NB_NODE-1);
  if (nb_extra_isets<=nb_conflict) {
    return NONE;
  }
  if (nb_extra_isets==nb_conflict+1 && 
      test_by_eliminate_failed_nodes()==NONE) {
    return NONE;
  }
    // printf("#%d#", ++NB1);
  return nb_conflict;
}

void printClique(int clique, int debut) {
  int i;
  printf("clique %d (%d nodes): ", clique, CLIQUE_STACK_fill_pointer-debut);
  for(i=debut; i<CLIQUE_STACK_fill_pointer; i++)
    printf("%d ", CLIQUE_STACK[i]);
  printf("\n");
}

int saved_current_MAX_CLQ_SIZE[tab_node_size];

void store_max_clique() {
  int i, node;
  MAXCLIQUE_STACK_fill_pointer=CLIQUE_STACK_fill_pointer;
  for(i=0; i<CLIQUE_STACK_fill_pointer; i++) 
    MAXCLIQUE_STACK[i]=CLIQUE_STACK[i];
  if (MAX_CLQ_SIZE != CLIQUE_STACK_fill_pointer)
    printf("erreur clique\n");
}

int init_MAX_CLQ_SIZE;

int solved_node_ub[tab_node_size*depth];
int candidate_ub[tab_node_size*depth];
int candidate_brother_debut[tab_node_size];
int candidate_brother_end[tab_node_size];

int re_number(int node) {
  int i, k, rest_clq_size, *adjacences, *adjacences1,
    *neibors, *saved_neibors, neibor, one_neibor;
  rest_clq_size=MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer;
  adjacences=matrice[node];
  for(i=0; i<rest_clq_size; i++) {
    neibors=ISETS[i]; one_neibor=NONE;
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) {
      if (adjacences[neibor]==TRUE) {
	if (one_neibor==NONE) {
	  one_neibor=neibor; saved_neibors=neibors;
	}
	else if (one_neibor!=NONE)
	  break;
      }
    }
    if (neibor==NONE) {
      adjacences1=matrice[one_neibor];
      for(k=i+1; k<rest_clq_size; k++) {
	neibors=ISETS[k];
	for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) {
	  if (adjacences1[neibor]==TRUE)
	    break;
	}
	if (neibor==NONE) {
	  ISETS[k][ISETS_SIZE[k]]=one_neibor;
	  ISETS_SIZE[k]++;
	  ISETS[k][ISETS_SIZE[k]]=NONE;
	  node_iset_no[one_neibor]=k;
	  *saved_neibors=node;
	  node_iset_no[node]=i;
	  return TRUE;
	}
      }
    }
  }
  return FALSE;
}

int my_addIntoIsetTomitaBis(int node) {
  int *neibors, neibor, i, *adjacences;
  adjacences=matrice[node];
  for(i=0; i<ISET_NB; i++) {
    neibors=ISETS[i];
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) {
      if (adjacences[neibor]==TRUE)
	break;
    }
    if (neibor==NONE) {
      if (i>MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer &&
	  i==ISET_NB-1)
	if (re_number(node)==TRUE)
	  return TRUE;
      ISETS[i][ISETS_SIZE[i]]=node;
      ISETS_SIZE[i]++;
      ISETS[i][ISETS_SIZE[i]]=NONE;
      node_iset_no[node]=i;
      return TRUE;
    }
  }
  if (i>MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer && re_number(node)==TRUE)
    return TRUE;
  else {
    ISETS_SIZE[ISET_NB]=1;
    ISETS[ISET_NB][0]=node;
    ISETS[ISET_NB][1]=NONE;
    node_iset_no[node]=ISET_NB;
    ISET_NB++;
    return FALSE;
  }  
}

int SHORT_ISET_STACK[tab_node_size];
int SHORT_ISET_STACK_fill_pointer=0;

int remove_node_from_iset(int node) {
  int iset, *noneibors, i, iset_size;

  iset=node_iset_no[node]; iset_size=ISETS_SIZE[iset];
  noneibors=ISETS[iset];
  for(i=0; i<iset_size; i++) {
    if (noneibors[i]==node) {
      iset_size=--ISETS_SIZE[iset];
      noneibors[i]=ISETS[iset][iset_size];
      ISETS[iset][iset_size]=NONE;
      if (iset_size==4) push(iset, SHORT_ISET_STACK);
      else if (iset_size==0) return NONE;
      return TRUE;
    }
  }
  printf("bizzar iset...");
  return FALSE;
}

int unitiset_filter(int iset) {
  int node, unitisetnode, *adjacences, i;
  unitisetnode=ISETS[iset][0];
  adjacences=matrice[unitisetnode];
  for(i=0;  i<FILTER_STACK_fill_pointer; i++) {
    node=FILTER_STACK[i];
    if (unitisetnode !=node && 
	filter_state[node]==ACTIVE && adjacences[node]!=TRUE) {
      filter_state[node]=NONE;
      if (remove_node_from_iset(node)==NONE)
	return NONE;
    }
  }
  return TRUE;
}

int biniset_filter(int iset) {
  int node, *adjacences1, *adjacences2, i;
  adjacences1=matrice[ISETS[iset][0]]; 
  adjacences2=matrice[ISETS[iset][1]];
  for(i=0;  i<FILTER_STACK_fill_pointer; i++) {
    node=FILTER_STACK[i];
    if (node_iset_no[node] != iset && filter_state[node]==ACTIVE && 
	adjacences1[node]!=TRUE && adjacences2[node]!=TRUE) {
      filter_state[node]=NONE;
      if (remove_node_from_iset(node)==NONE)
	return NONE;
    }
  }
  return TRUE;
}

int triset_filter(int iset) {
  int i, node, *adjacences1, *adjacences2, *adjacences3;
  adjacences1=matrice[ISETS[iset][0]]; 
  adjacences2=matrice[ISETS[iset][1]];
  adjacences3=matrice[ISETS[iset][2]];
  for(i=0;  i<FILTER_STACK_fill_pointer; i++) {
    node=FILTER_STACK[i];
    if (node_iset_no[node] != iset && filter_state[node]==ACTIVE && 
	adjacences1[node]!=TRUE && adjacences2[node]!=TRUE &&
	adjacences3[node]!=TRUE) {
      filter_state[node]=NONE;
      if (remove_node_from_iset(node)==NONE)
	return NONE;
    }
  }
  return TRUE;
}

int fouriset_filter(int iset) {
  int i, node, *adjacences1, *adjacences2, *adjacences3,
    *adjacences4;
  adjacences1=matrice[ISETS[iset][0]]; 
  adjacences2=matrice[ISETS[iset][1]];
  adjacences3=matrice[ISETS[iset][2]];
  adjacences4=matrice[ISETS[iset][3]];
  for(i=0;  i<FILTER_STACK_fill_pointer; i++) {
    node=FILTER_STACK[i];
    if (node_iset_no[node] != iset && filter_state[node]==ACTIVE && 
	adjacences1[node]!=TRUE && adjacences2[node]!=TRUE &&
	adjacences3[node]!=TRUE && adjacences4[node]!=TRUE) {
      filter_state[node]=NONE;
      if (remove_node_from_iset(node)==NONE)
	return NONE;
    }
  }
  return TRUE;
}

int filter() {
  int iset, i;
  SHORT_ISET_STACK_fill_pointer=0;
  for(i=0; i<ISET_NB; i++)
    if (ISETS_SIZE[i]<=4)
      push(i, SHORT_ISET_STACK);
  for(i=0; i<SHORT_ISET_STACK_fill_pointer; i++) {
    iset=SHORT_ISET_STACK[i];
    if (ISETS_SIZE[iset]==1 && unitiset_filter(iset)==NONE)
      return NONE;
    else if (ISETS_SIZE[iset]==2 && biniset_filter(iset)==NONE)
      return NONE;
    else if (ISETS_SIZE[iset]==3 && triset_filter(iset)==NONE)
      return NONE;
    else if (ISETS_SIZE[iset]==4 && fouriset_filter(iset)==NONE)
      return NONE;
  }
  return TRUE;
}

void eliminate_incompatibe_nodes() {
  int node, i, j;
  for(i=0; i<FILTER_STACK_fill_pointer; i++) {
    if (filter_state[FILTER_STACK[i]]==NONE) {
      break;
    }
  }
  for(j=i+1; j<FILTER_STACK_fill_pointer; j++) {
    node=FILTER_STACK[j];
    if (filter_state[node]==ACTIVE)
      FILTER_STACK[i++]=node;
  }
  FILTER_STACK_fill_pointer=i;
}

int my_get_upper_bound(int candidate, int debut) {
  int i, node;

  NB_ISETS_TEST++; 
  ISET_NB=0; 
  for(i=0; i<FILTER_STACK_fill_pointer; i++) {
    node=FILTER_STACK[i];
    my_addIntoIsetTomitaBis(node);
    filter_state[node]=ACTIVE;
  }
  if (ISET_NB==MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer) {
    if (filter()==NONE) 
      return ISET_NB-1;
    //    else eliminate_incompatibe_nodes();
  }
  return ISET_NB;
}

int my_maxsatz(int node, int rest_clq_size) {
  int i, *neibors, neibor, nb_conflict;
  ISETS_SIZE[ISET_NB]=1;
  ISETS[ISET_NB][0]=node;
  ISETS[ISET_NB][1]=NONE;
  node_iset_no[node]=ISET_NB;
  ISET_NB++;
  for(i=0; i<FILTER_STACK_fill_pointer; i++) 
    if (filter_state[FILTER_STACK[i]]==ACTIVE)
      node_state[FILTER_STACK[i]]=ACTIVE;
  node_state[node]=ACTIVE;
  nb_conflict=maxsatz(rest_clq_size);
  reset_enlarged_isets();
  for(i=0; i<FILTER_STACK_fill_pointer; i++) 
    node_state[FILTER_STACK[i]]=PASSIVE;
  node_state[node]=PASSIVE;
  ISET_NB--;
  if (nb_conflict==NONE) return rest_clq_size;
  else return ISET_NB-nb_conflict+1;
}

int estimate_ub(int node, int old_ub, int debut) {
  int i, max=0, reduced_nb, rest_clq_size, ub;
  FILTER_STACK_fill_pointer=0;
  for(i=debut; i<SOLVED_NODE_STACK_fill_pointer; i++) {
    if (matrice[node][SOLVED_NODE_STACK[i]]==TRUE) {
      if (max<solved_node_ub[i])
	max=solved_node_ub[i];
      push(SOLVED_NODE_STACK[i], FILTER_STACK);
    }
  }
  max++; // add node itself into the clique
  if (old_ub<max)
    max=old_ub;
  rest_clq_size=MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer;
  if (max>rest_clq_size) {
    ub=my_get_upper_bound(node, debut);
    if (max>ub+1)
      max=ub+1;
    if (max>rest_clq_size) {
      reduced_nb=my_maxsatz(node, rest_clq_size);
      if (max>reduced_nb)
	max=reduced_nb;
    }
  }
  return max;
}

int estimate_ub2(int node, int ub, int debut) {
  int i;
  for(i=debut-1; i>=0; i--) {
    if (SOLVED_NODE_STACK[i]==node) {
      if (solved_node_ub[i]<ub)
	return solved_node_ub[i];
      else return ub;
    }
  }
  return ub;
}

int backtracking_for_maxclique(int *debut) {
  int node;
  NB_BACK_CLIQUE++;
  node=pop(CLIQUE_STACK);
  while (CLIQUE_CANDIDATE_STACK_fill_pointer>0 &&
	 CLIQUE_CANDIDATE_STACK[CLIQUE_CANDIDATE_STACK_fill_pointer-1]==NONE) {
    pop(CLIQUE_CANDIDATE_STACK);
    node=pop(CLIQUE_STACK);
  }
  SOLVED_NODE_STACK_fill_pointer=candidate_brother_end[node];
  *debut=candidate_brother_debut[node];
  solved_node_ub[SOLVED_NODE_STACK_fill_pointer]=MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer;
  push(node, SOLVED_NODE_STACK);
  return candidate_brother_debut[node];
}

int expand_clq_from_node1(int node) {
  int i, nb, candidate, debut, *neibors;
  debut=candidate_brother_debut[node]; neibors=matrice[node]; nb=0;
  for(i=SOLVED_NODE_STACK_fill_pointer-1; i>=debut; i--) {
    candidate=SOLVED_NODE_STACK[i];
    if (neibors[candidate]==TRUE && filter_state[candidate]==ACTIVE) {
      nb++;
      candidate_ub[CLIQUE_CANDIDATE_STACK_fill_pointer]=solved_node_ub[i];
      push(candidate, CLIQUE_CANDIDATE_STACK);
    }
  }
  return nb;
}

int estimate_ub3(int node, int old_ub, int debut) {
  int i, max=0, max_iset_nb=0, solved_node;
  for(i=debut; i<SOLVED_NODE_STACK_fill_pointer; i++) {
    solved_node=SOLVED_NODE_STACK[i];
    if (max_iset_nb<node_iset_no[solved_node])
      max_iset_nb=node_iset_no[solved_node];
    if (matrice[node][SOLVED_NODE_STACK[i]]==TRUE) {
      if (max<solved_node_ub[i])
	max=solved_node_ub[i];
    }
  }
  max++; // add node itself into the clique
  if (old_ub<max)
    max=old_ub;
  if (max_iset_nb<node_iset_no[node])
    max_iset_nb=node_iset_no[node];
  if (max_iset_nb+1<max)  // the number of iset begins from 0
    max=max_iset_nb+1;
  return max;
}

#define NO_NEW_ISET -7
#define NEW_ISET -77

int min2(int a, int b) {
  if (a<b)
    return a;
  else return b;
}

void insert_node_into_conflict_iset(int node, int iset) {
  int *neibors, neibor, next_neibor;
  neibors=ISETS[iset];
  for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors))
    if (neibor>NB_NODE)
      break;
  *neibors=node; node_iset_no[node]=iset;
  for(next_neibor=*(++neibors); next_neibor!=NONE; 
      next_neibor=*(++neibors)) {
    *neibors=neibor, neibor=next_neibor;
  }
  *neibors=neibor, *(neibors+1)=NONE;
  ISETS_SIZE[iset]++;
}

void insert_node_into_iset(int node, int i) {
  ISETS[i][ISETS_SIZE[i]]=node;
  ISETS_SIZE[i]++;
  ISETS[i][ISETS_SIZE[i]]=NONE;
  node_iset_no[node]=i;
}

void open_new_iset_for_node(int node) {
  ISETS_SIZE[ISET_NB]=1;
  ISETS[ISET_NB][0]=node;
  ISETS[ISET_NB][1]=NONE;
  node_iset_no[node]=ISET_NB;
  iset_state[ISET_NB]=ACTIVE;
  ISET_NB++;
}

// the function returns no_new_iset if node is inserted into
// a non-relaxed iset, new_iset if a new iset is opened for node
// (the new iset is unit now), or the last relaxing variable
// in the relaxed iset in which node is inserted.
int my_inc_addIntoIsetTomitaBis(int node) {
  int *neibors, neibor, i, *adjacences;
  // relaxed,  last_relaxing_var=NONE, relaxed_iset;
  adjacences=matrice[node];
  for(i=0; i<ISET_NB; i++) {
    neibors=ISETS[i]; // relaxed=FALSE;
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) {
      if (neibor<=NB_NODE) {
	if (adjacences[neibor]==TRUE)
	  break;
      }
      else {
	break;
/* 	relaxed=TRUE; */
/* 	if (last_relaxing_var==NONE) { */
/* 	  for(last_relaxing_var=neibor, neibor=*(++neibors);   */
/* 	      neibor!=NONE; last_relaxing_var=neibor, neibor=*(++neibors)); */
/* 	  relaxed_iset=i; */
/* 	} */
      }
    }
    if (neibor==NONE) { // && relaxed==FALSE) {
      insert_node_into_iset(node, i); return NO_NEW_ISET;
    }
  }
  open_new_iset_for_node(node); return NEW_ISET;
/*   if (last_relaxing_var==NONE) { */
/*     open_new_iset_for_node(node); return NEW_ISET; */
/*   } */
/*   else { */
/*     open_new_iset_for_node(node); return NEW_ISET; */
/*     insert_node_into_conflict_iset(node, relaxed_iset); */
/*     return last_relaxing_var; */
}

void init_inc_maxsatz() {
  int i;

  UNITISET_STACK_fill_pointer=0;
  for(i=0; i<ISET_NB; i++) {
    iset_used[i]=FALSE;
    iset_state[i]=ACTIVE;
    if (ISETS_SIZE[i]==1)
      push(i, UNITISET_STACK);
  }
  ADDED_NB_NODE=NB_NODE+1;
  CONFLICT_ISET_STACK_fill_pointer=0;
  ENLARGED_ISET_STACK_fill_pointer=0;
}

int inc_test_node(int node, int iset, int further) {
  int my_iset, saved_unitiset_stack_fill_pointer, 
    saved_node_stack_fill_pointer, saved_passive_iset_stack_fill_pointer,
    saved_reduced_iset_stack_fill_pointer;

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  saved_node_stack_fill_pointer=NODE_STACK_fill_pointer;
  saved_passive_iset_stack_fill_pointer=PASSIVE_ISET_STACK_fill_pointer;
  saved_reduced_iset_stack_fill_pointer=REDUCED_ISET_STACK_fill_pointer;
  MY_UNITISET_STACK_fill_pointer=0;
  //  if (ISETS_SIZE[iset]>1 && (my_iset=unitIsetProcess())!=NO_CONFLICT)
  //  printf("bizzar.....\n");
  if ((my_iset=fix_node_for_iset(node, iset)) != NO_CONFLICT ||
      (my_iset=my_unitIsetProcess()) !=  NO_CONFLICT ||
      (further &&
       (my_iset=further_test_node(saved_reduced_iset_stack_fill_pointer)) 
       != NO_CONFLICT)) {
    lookback_for_maxsatz(my_iset);
    reset_context_for_maxsatz(saved_node_stack_fill_pointer,
			      saved_passive_iset_stack_fill_pointer,
			      saved_reduced_iset_stack_fill_pointer,
			      saved_unitiset_stack_fill_pointer);
    return my_iset;
  }
  else {
    reset_context_for_maxsatz(saved_node_stack_fill_pointer,
				 saved_passive_iset_stack_fill_pointer,
				 saved_reduced_iset_stack_fill_pointer,
				 saved_unitiset_stack_fill_pointer);
    return NO_CONFLICT;
  }
}

int inc_maxsatz_lookahead_by_fl(int saved_unitiset_stack_fill_pointer) {
  int iset, *nodes, node, no_conflict, test_flag, k, saved_size, ctr, 
    i, one_conflict=FALSE, result;
  for(k=2; k<=FL_TEST_LENGTH; k++) 
    //for(i=0; i<=ISET_NB-nb_extra_isets+nb_conflict; i++) {
   for(iset=ISET_NB-1; iset>=0; iset--) {
    if (iset_state[iset]==ACTIVE && ISETS_SIZE[iset]==k) {
      nodes=ISETS[iset]; test_flag=TRUE;
      for(node=*nodes; node!=NONE; node=*(++nodes)) {
	if (node>NB_NODE && node_state[node]==ACTIVE) {
	  test_flag=FALSE; break;
	}
      }
      if (test_flag==TRUE) {
	nodes=ISETS[iset]; no_conflict=FALSE; ctr=k+1;
	INVOLVED_ISET_STACK_fill_pointer=0;
	for(node=*nodes; node!=NONE; node=*(++nodes))
	  if (node_state[node]==ACTIVE) {
	    if (inc_test_node(node, iset, (*(nodes+1)==NONE))==NO_CONFLICT) {
	      no_conflict=TRUE; break;
	    }
	    else store_involved_isets();
	  }
	if (no_conflict==FALSE) { 
	  reset_context_for_maxsatz(0, 0, 0,
			    saved_unitiset_stack_fill_pointer);
	  saved_size=ISETS_SIZE[iset];
	  enlarge_stored_involved_isets();
	  if (saved_size+1 != ISETS_SIZE[iset])
	    printf("erreur iset involved %d %d %d...", 
		   saved_size,  ISETS_SIZE[iset], iset);
	  one_conflict=TRUE;
	}
	for(i=0; i<INVOLVED_ISET_STACK_fill_pointer; i++)
	  iset_used[INVOLVED_ISET_STACK[i]]=FALSE;
	INVOLVED_ISET_STACK_fill_pointer=0;
	if (one_conflict==TRUE) 
	  return NONE;
      }
    }
  }
  return TRUE;
}

int my_inc_maxsatz(int node) {
  int iset, saved_unitiset_stack_fill_pointer, result, my_iset;

  saved_unitiset_stack_fill_pointer=UNITISET_STACK_fill_pointer;
  NODE_STACK_fill_pointer=0; 
  PASSIVE_ISET_STACK_fill_pointer=0;
  REDUCED_ISET_STACK_fill_pointer=0;
  my_iset=UNITISET_STACK[0];
  UNITISET_STACK[0]=UNITISET_STACK[UNITISET_STACK_fill_pointer-1];
  UNITISET_STACK[UNITISET_STACK_fill_pointer-1]=my_iset;
  MY_UNITISET_STACK_fill_pointer=0;
  
  if ((iset=unitIsetProcess())!=NO_CONFLICT) {
    lookback_for_maxsatz(iset);
    reset_context_for_maxsatz(0, 0, 0,
			      saved_unitiset_stack_fill_pointer);
    enlarge_involved_iset();
    return NONE;
  }
  result=inc_maxsatz_lookahead_by_fl(saved_unitiset_stack_fill_pointer);
  reset_context_for_maxsatz(0, 0, 0,
			    saved_unitiset_stack_fill_pointer);
  if (result==NONE)
    return NONE;
  return FALSE;
}

int NB_INC_SUCCESS=0;

int solve_first_nodes(int debut) {
  int i, node, rest_clq_size, old_iset_nb, result, init_maxsatz=FALSE;
  ISET_NB=0; rest_clq_size=MAX_CLQ_SIZE-CLIQUE_STACK_fill_pointer;
  FILTER_STACK_fill_pointer=0;
  while (CLIQUE_CANDIDATE_STACK_fill_pointer>0) {
    node=pop(CLIQUE_CANDIDATE_STACK);
    if (node==NONE) {
      NB_INC_SUCCESS++;
      //   printf("jhgjhg "); NB_INC_SUCCESS++;
      break;
    }
    else {   node_state[node]=ACTIVE;
      push(node, FILTER_STACK);
      result=my_inc_addIntoIsetTomitaBis(node);
      if (ISET_NB<=rest_clq_size || result==NO_NEW_ISET) {
	solved_node_ub[SOLVED_NODE_STACK_fill_pointer]=
	  estimate_ub3(node,  min2(rest_clq_size, 
			  candidate_ub[CLIQUE_CANDIDATE_STACK_fill_pointer]),
		       debut);
	push(node, SOLVED_NODE_STACK);
      }
      else if (result==NEW_ISET) {
	if (init_maxsatz==FALSE) {
	  init_maxsatz=TRUE; init_inc_maxsatz();
	}   
	push(ISET_NB-1, UNITISET_STACK); //the new iset is ISET_NB-1 and is unit
	if (my_inc_maxsatz(node)==NONE) {
	  solved_node_ub[SOLVED_NODE_STACK_fill_pointer]=
	    estimate_ub3(node, min2(rest_clq_size, 
			  candidate_ub[CLIQUE_CANDIDATE_STACK_fill_pointer]), 
			 debut);
	  push(node, SOLVED_NODE_STACK);
	}
	else break;
      }
      else {
 // result is an added node to relax an iset in a subset of conflicting isets
/* 	if (my_inc_maxsatz2(node, result)==NONE) { */
/* 	  solved_node_ub[SOLVED_NODE_STACK_fill_pointer]= */
/* 	    estimate_ub3(node,   */
/* 			 min2(rest_clq_size,  */
/* 			    candidate_ub[CLIQUE_CANDIDATE_STACK_fill_pointer]),  */
/* 			 debut); */
/* 	  push(node, SOLVED_NODE_STACK); */
/* 	} */
/* 	else */ break;
      }
    }
  }
  for(i=debut; i<SOLVED_NODE_STACK_fill_pointer; i++)
      node_state[SOLVED_NODE_STACK[i]]=PASSIVE;
  if (node!=NONE) node_state[node]=PASSIVE;
  reset_enlarged_isets();
  push(node, CLIQUE_CANDIDATE_STACK);
  return TRUE;
}


void printMaxClique() {
  int i;
  printf("Max Clique Size: %d\n", MAX_CLQ_SIZE);
  printf("Clique Output:");
  for(i=0; i<MAXCLIQUE_STACK_fill_pointer; i++)
    printf(" %d", MAXCLIQUE_STACK[i]);
  printf("\n");
}

int maxClique(int cutoff) {
  int i, node, candidate, highdensity, res, debut=0, extra_iset_nb, iset_nb, ub;
  int min, j, *neibors, neibor;

  CLIQUE_STACK_fill_pointer=0; 
  SOLVED_NODE_STACK_fill_pointer=0; NB_BACK_CLIQUE=0;
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    candidate_ub[i]=NB_NODE;
    node_state[CLIQUE_CANDIDATE_STACK[i]]=PASSIVE;
  }
  while (CLIQUE_CANDIDATE_STACK_fill_pointer>0) {
    candidate=pop(CLIQUE_CANDIDATE_STACK);
    if (candidate==NONE) {
      backtracking_for_maxclique(&debut);
      if (NB_BACK_CLIQUE>cutoff)
	break;
    }
    else {
      ub=estimate_ub(candidate,
		     candidate_ub[CLIQUE_CANDIDATE_STACK_fill_pointer], debut);
      if (CLIQUE_STACK_fill_pointer+ub<=MAX_CLQ_SIZE) {
	solved_node_ub[SOLVED_NODE_STACK_fill_pointer]=ub;
	push(candidate, SOLVED_NODE_STACK);
      }
      else {
	saved_candidate_stack[candidate]=CLIQUE_CANDIDATE_STACK_fill_pointer;
	push(NONE, CLIQUE_CANDIDATE_STACK);
	push(candidate, CLIQUE_STACK);
	candidate_brother_debut[candidate]=debut;
	candidate_brother_end[candidate]=SOLVED_NODE_STACK_fill_pointer;
	debut=SOLVED_NODE_STACK_fill_pointer;
	res=expand_clq_from_node1(candidate);
	if (res==0) {
	  if (CLIQUE_STACK_fill_pointer>MAX_CLQ_SIZE) {
	    MAX_CLQ_SIZE=CLIQUE_STACK_fill_pointer;
	    store_max_clique();
	    printf("Current MaxClique Size=%d %d\n", 
		   MAX_CLQ_SIZE, NB_BACK_CLIQUE);
            printMaxClique();
	    node=CLIQUE_STACK[0];
	    CLIQUE_STACK_fill_pointer=0;
	    CLIQUE_CANDIDATE_STACK_fill_pointer=saved_candidate_stack[node];
	    debut=0;
	    SOLVED_NODE_STACK_fill_pointer=candidate_brother_end[node];
	    solved_node_ub[SOLVED_NODE_STACK_fill_pointer]=MAX_CLQ_SIZE;
	    push(node, SOLVED_NODE_STACK);
	  }
	  else 
	    printf("bizzar....MAX_CLQ_SIZE=%d Current clique size=%d\n", 
		   MAX_CLQ_SIZE, CLIQUE_STACK_fill_pointer);
/* 	  pop(CLIQUE_CANDIDATE_STACK); */
/* 	  backtracking_for_maxclique(&debut); */
	}
	else solve_first_nodes(debut);
      }
    }
  }
  return MAX_CLQ_SIZE;
}

int NB_CLIQUE;

static int degreecmp_static_inc(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degreeNeibors1, degreeNeibors2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=node_nb_neibors[*node1]; degree2=node_nb_neibors[*node2];
  degreeNeibors1=degree_neibors[*node1];
  degreeNeibors2=degree_neibors[*node2];
  if (degree1>degree2)
    return 1;
  else if (degree1==degree2) {
    if (degreeNeibors1>degreeNeibors2)
      return 1;
    else if (degreeNeibors1==degreeNeibors2)
      return 0;
    else return -1;
  }
  else return -1;
}

void sort_by_degree() {
  int min, i, node, neibor, *neibors, j;
  CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for(i=1; i<=NB_NODE; i++) {
    min=NB_NODE;
    for(j=1; j<=NB_NODE; j++) {
      if (node_state[j]==ACTIVE && node_nb_neibors[j]<min) { 
	min=node_nb_neibors[j]; node=j;
      }
    }
    push(node, CLIQUE_CANDIDATE_STACK);
    node_state[node]=PASSIVE;
    neibors=node_neibors[node];
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) 
      if (node_state[neibor]==ACTIVE) {
	node_nb_neibors[neibor]--;
      }
  }
}

int init_for_maxclique_again() {
    
  int node, i;
  MAX_CLQ_SIZE=iClique(0);
  store_max_clique();
  
  init_MAX_CLQ_SIZE=MAX_CLQ_SIZE;
  for(i=0; i<CLIQUE_STACK_fill_pointer; i++)
    node_state[CLIQUE_STACK[i]]=ACTIVE;
  CLIQUE_STACK_fill_pointer=0;

  CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for(node=1; node<=NB_NODE; node++) {
   if (node_state[node]==ACTIVE) {
      push(node, CLIQUE_CANDIDATE_STACK);
    }
    nodeMark[node]=0;
    filter_state[node]=PASSIVE;
    iset_state[node]=ACTIVE;
    iset_involved[node]=FALSE;
    //  conflict_iset[node]=FALSE;
    node_reason[node]=NO_REASON;
    node_tested_state[node]=FALSE;
  }
  NONEIBOR_STACK_fill_pointer=0;
  CLIQUE_STACK_fill_pointer=0;
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    candidate_potentiel[i]=i+1;
    //   candidate_nb_neibors[i]=node_nb_neibors[CLIQUE_CANDIDATE_STACK[i]];
  }
  //compute_candidate_nb_neibors(0);
  //compute_degree_neibors(0);
  for(node=1; node<=NB_NODE; node++) {
    nb_neibors[node]=static_nb_neibors[node];
    node_nb_neibors[node]=static_nb_neibors[node];
  }
  qsort(CLIQUE_CANDIDATE_STACK, 
	CLIQUE_CANDIDATE_STACK_fill_pointer, 
	sizeof(int), degreecmp_static_inc);
  sort_by_degree();
}

void my_sort_by_degree_dec() {
  int max, i, node, neibor, *neibors, j, node1, index;
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i]; max=nb_neibors[node]; index=i;
    for(j=i+1; j<CLIQUE_CANDIDATE_STACK_fill_pointer; j++) {
      node1=CLIQUE_CANDIDATE_STACK[j];
      if (nb_neibors[node1]>max) { 
	node=node1; max=nb_neibors[node1]; index=j;
      }
    }
    node1=CLIQUE_CANDIDATE_STACK[i];
    CLIQUE_CANDIDATE_STACK[i]=node;
    CLIQUE_CANDIDATE_STACK[index]=node1;
    neibors=node_neibors[node];
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) 
      nb_neibors[neibor]--;
  }
/*   for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) { */
/*     node=CLIQUE_CANDIDATE_STACK[i]; */
  for(node=1; node<=NB_NODE; node++) {
    nb_neibors[node]=static_nb_neibors[node];
  }
}

void my_sort_by_degree_inc() {
  int min, i, node, neibor, *neibors, j, node1, index;
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    node=CLIQUE_CANDIDATE_STACK[i]; min=nb_neibors[node]; index=i;
    for(j=i+1; j<CLIQUE_CANDIDATE_STACK_fill_pointer; j++) {
      node1=CLIQUE_CANDIDATE_STACK[j];
      if (nb_neibors[node1]<min) { 
	node=node1; min=nb_neibors[node1]; index=j;
      }
    }
    node1=CLIQUE_CANDIDATE_STACK[i];
    CLIQUE_CANDIDATE_STACK[i]=node;
    CLIQUE_CANDIDATE_STACK[index]=node1;
    neibors=node_neibors[node];
    for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) 
      nb_neibors[neibor]--;
  }
/*   for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) { */
/*     node=CLIQUE_CANDIDATE_STACK[i]; */
  for(node=1; node<=NB_NODE; node++) {
    nb_neibors[node]=static_nb_neibors[node];
  }
}

double ISETS_SumDegree[tab_node_size];
static int degreecmp_dec(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degreeNeibors1, degreeNeibors2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=node_nb_neibors[*node1]; degree2=node_nb_neibors[*node2];
  if (degree1>degree2)
    return -1;
  else if (degree1==degree2) {
      return 0;
  }
  else return 1;
}

int iset_smaller_than(int iset1, int iset2) {
  int *neibors1, neibor1, *neibors2, neibor2;
/*   if (ISETS_SumDegree[iset1]<ISETS_SumDegree[iset2]) */
/*     return TRUE; */
/*   else if (ISETS_SumDegree[iset1]>ISETS_SumDegree[iset2]) */
/*     return FALSE; */
  if (ISETS_SIZE[iset1]<ISETS_SIZE[iset2])
    return TRUE;
  else if (ISETS_SIZE[iset1]>ISETS_SIZE[iset2])
    return FALSE;
  neibors1=ISETS[iset1];   neibors2=ISETS[iset2];
  for(neibor1=*neibors1, neibor2=*neibors2; neibor1 !=NONE && neibor2 !=NONE;
      neibor1=*(++neibors1),neibor2=*(++neibors2)) {
    if (node_nb_neibors[neibor1]<node_nb_neibors[neibor2])
      return TRUE;
    else if (node_nb_neibors[neibor1]>node_nb_neibors[neibor2])
      return FALSE;
  }
  if (neibor1==NONE && neibor2==NONE)
    return FALSE;
  else if (neibor1==NONE)
    return TRUE;
  else return FALSE;
}


int sort_isets_and_push_nodes() {
  int i, j, *neibors, neibor, iset1, iset2, node, index;
  ISET_STACK_fill_pointer=0; CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for(i=0; i<ISET_NB; i++) {
    push(i, ISET_STACK);
    neibors=ISETS[i]; ISETS_SumDegree[i]=1.0;
    for(neibor=*neibors; neibor !=NONE; neibor=*(++neibors))
      ISETS_SumDegree[i]=ISETS_SumDegree[i]*node_nb_neibors[neibor];
  }
  for(i=0; i<ISET_STACK_fill_pointer; i++) {
    iset1=ISET_STACK[i]; index=i;
    qsort(ISETS[iset1], ISETS_SIZE[iset1], sizeof(int), degreecmp_dec);
    for (j=i+1; j<ISET_STACK_fill_pointer; j++) {
      iset2=ISET_STACK[j];
      qsort(ISETS[iset2], ISETS_SIZE[iset2], sizeof(int), degreecmp_dec);
      if (iset_smaller_than(iset2, iset1)==TRUE) {
	iset1=iset2; index=j;
      }
    }
    iset2=ISET_STACK[i]; ISET_STACK[i]=iset1, ISET_STACK[index]=iset2;
    for(j=ISETS_SIZE[iset1]-1; j>=0; j--) {
      node=ISETS[iset1][j];
      push(node, CLIQUE_CANDIDATE_STACK);
/*       neibors=node_neibors[node]; */
/*       for(neibor=*neibors; neibor!=NONE; neibor=*(++neibors)) { */
/* 	node_nb_neibors[neibor]--; */
/*       } */
    }
  }
}


int sort_by_maxIsets() {
  int node, i, nb_isets=0, ordered_nodes[tab_node_size], nb_node;
  complement_graph();
  CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for(node=1; node<=NB_NODE; node++)
    if (node_state[node] !=NONE)
      push(node, CLIQUE_CANDIDATE_STACK);
  if (NB_NODE<1000)
    my_sort_by_degree_dec();
  else
    my_sort_by_degree_inc();
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++)
    ordered_nodes[i]=CLIQUE_CANDIDATE_STACK[i];
  nb_node=CLIQUE_CANDIDATE_STACK_fill_pointer;
  int ans = 0;
  while (CLIQUE_CANDIDATE_STACK_fill_pointer>0) {
    //qsort(CLIQUE_CANDIDATE_STACK, CLIQUE_CANDIDATE_STACK_fill_pointer, 
    //	  sizeof(int), degreecmp_dec);
    MAX_CLQ_SIZE=0;
    maxClique(50000);
    for(i=0; i<MAXCLIQUE_STACK_fill_pointer; i++) {
      node=MAXCLIQUE_STACK[i]; node_state[node]=NONE;
      push(node, ISET_STACK);
    }
    push(NONE, ISET_STACK); nb_isets++;
    CLIQUE_CANDIDATE_STACK_fill_pointer=0;
    for(i=0; i<nb_node; i++) {
      if (node_state[ordered_nodes[i]] !=NONE)
	push(ordered_nodes[i], CLIQUE_CANDIDATE_STACK);
    }
    if (MAXCLIQUE_STACK_fill_pointer==1) {
        ans++;
    }
    printf("****%d iset size: %d, nb back: %d #remaining nodes: %d\n", nb_isets, 
	   MAXCLIQUE_STACK_fill_pointer, NB_BACK_CLIQUE, 
	   CLIQUE_CANDIDATE_STACK_fill_pointer);
  }
  printf("init nb isets %d\n\n", nb_isets);
  complement_graph(); MAX_CLQ_SIZE=init_MAX_CLQ_SIZE;
  ISET_NB=0; ISETS_SIZE[ISET_NB]=0;
  // for(i=ISET_STACK_fill_pointer-1; i>=0; i--)
  for(i=0; i<ISET_STACK_fill_pointer; i++) {
    if (ISET_STACK[i]!=NONE)
      ISETS[ISET_NB][ISETS_SIZE[ISET_NB]++]=ISET_STACK[i];
    else {
      ISETS[ISET_NB][ISETS_SIZE[ISET_NB]]=NONE;
      ISET_NB++; ISETS_SIZE[ISET_NB]=0;
    }
  }
  sort_isets_and_push_nodes();
  if (ans>1) return FALSE;
  else return TRUE;
}



int init_for_maxclique(int ordering) {
  DENSITY=NB_EDGE*100*2/(NB_NODE*(NB_NODE-1));
//  RATIO=1.0*NB_NODE*(NB_NODE-1)/(2*NB_EDGE);
  int old = ordering;
  if (ordering ==1) 
    printf("Using the degeneracy ordering...\n");
  else if (ordering ==2)
    printf("Using the MaxIndSet ordering...\n");
  else if (ordering == -1) {
        if (DENSITY <= 70) {
	  printf("Using the degeneracy ordering...\n");
	  ordering = 1;
        } else {
	  printf("Investigating the MaxIndSet ordering...\n");
	  ordering = 2;
        }
  }
    
  int node, i;
  if (ordering==1) {
    MAX_CLQ_SIZE=iClique(0);
    store_max_clique();
  }
  else MAX_CLQ_SIZE=0;
  init_MAX_CLQ_SIZE=MAX_CLQ_SIZE;
  for(i=0; i<CLIQUE_STACK_fill_pointer; i++)
    node_state[CLIQUE_STACK[i]]=ACTIVE;
  CLIQUE_STACK_fill_pointer=0;

  CLIQUE_CANDIDATE_STACK_fill_pointer=0;
  for(node=1; node<=NB_NODE; node++) {
   if (node_state[node]==ACTIVE) {
      push(node, CLIQUE_CANDIDATE_STACK);
    }
    nodeMark[node]=0;
    filter_state[node]=PASSIVE;
    iset_state[node]=ACTIVE;
    iset_involved[node]=FALSE;
    //  conflict_iset[node]=FALSE;
    node_reason[node]=NO_REASON;
    node_tested_state[node]=FALSE;
  }
  NONEIBOR_STACK_fill_pointer=0;
  CLIQUE_STACK_fill_pointer=0;
  for(i=0; i<CLIQUE_CANDIDATE_STACK_fill_pointer; i++) {
    candidate_potentiel[i]=i+1;
    //   candidate_nb_neibors[i]=node_nb_neibors[CLIQUE_CANDIDATE_STACK[i]];
  }
  compute_candidate_nb_neibors(0);
  compute_degree_neibors(0);

  qsort(CLIQUE_CANDIDATE_STACK, 
	CLIQUE_CANDIDATE_STACK_fill_pointer, 
	sizeof(int), degreecmp_static_inc);
  if (ordering==1)
    sort_by_degree();
  else {
    int ans = sort_by_maxIsets();
    if (ans == FALSE && old == -1) {
        int node = 1;
        for (; node <= NB_NODE; ++node) {
            node_state[node] = ACTIVE;
        }
        printf("the MaxIndSet ordering disable, using the degeneracy ordering...\n");
        init_for_maxclique_again();
        return 1;
    }
    else if (old ==-1) printf("Using the MaxIndSet ordering...\n");
  }
  return ordering;
}


static int degreecmp_inc(const void *pnode1, const void *pnode2) {
  int *node1, *node2, degree1, degree2, degreeNeibors1, degreeNeibors2;
  //  node1=(*pnode1); node2=(*pnode2);
  node1=(int *) pnode1; node2=(int *) pnode2;
  degree1=node_nb_neibors[*node1]; degree2=node_nb_neibors[*node2];
  if (degree1>degree2)
    return 1;
  else if (degree1==degree2) {
      return 0;
  }
  else return -1;
}

//int ISET_STACK[tab_node_size], ISET_STACK_fill_pointer=0;
int iset_degree[tab_node_size];
int iset_degree2[tab_node_size];

static int isetDegreeCMP(const void *p_iset1, const void *p_iset2) {
  int *iset1, *iset2;
  iset1=(int *)p_iset1;   iset2=(int *)p_iset2; 
  if ((iset_degree[*iset1]<iset_degree[*iset2])
      || (iset_degree[*iset1]==iset_degree[*iset2] && 
	  iset_degree2[*iset1]<iset_degree2[*iset2]))
    return -1;
  else if (iset_degree[*iset1]==iset_degree[*iset2]
	   &&  iset_degree2[*iset1]==iset_degree2[*iset2])
    return 0;
  else return 1;
}

//double ISETS_SumDegree[tab_node_size];

main(int argc, char *argv[]) {
  int i, result, ordering;
  //  unsigned long long begintime, endtime, mess;
  //  struct tms *a_tms;
  FILE *fp_time;
  struct rusage starttime, endtime;
  
/*   a_tms = ( struct tms *) malloc( sizeof (struct tms)); */
/*   mess=times(a_tms); begintime = a_tms->tms_utime; */

  getrusage(RUSAGE_SELF, &starttime );

  if (argc>2) ordering=atoi(argv[2]);
  else ordering=-1;
  
  // begintime=clock();
  switch (build_simple_graph_instance(argv[1])) {
  case FALSE: printf("Input file error\n"); return FALSE;
  case TRUE:
    build_complement_graph();
    init_for_maxclique(ordering);
    printf("instance information: #node=%d, #edge=%d density= %5.4f\n", 
	   NB_NODE, NB_EDGE, ((float)NB_EDGE*2)/(NB_NODE*(NB_NODE-1)));
    printf("initial clique size: %d\n\n", MAX_CLQ_SIZE);
    maxClique(2000000000);
    printMaxClique();
    break;
  }

  //  mess=times(a_tms); endtime = a_tms->tms_utime;
  // endtime=clock();

    getrusage( RUSAGE_SELF, &endtime );

  printf ("Program terminated in %d seconds.\n",
	  ((int)(endtime.ru_utime.tv_sec - starttime.ru_utime.tv_sec)));

  printf("NewIncMaxCLQ %s %d %d %5.4f %d %d %d %d %d\n", 
	 argv[1], NB_NODE, NB_EDGE, 
	 ((float)NB_EDGE*2)/(NB_NODE*(NB_NODE-1)),
	 ((int)(endtime.ru_utime.tv_sec - starttime.ru_utime.tv_sec)), 
	 NB_BACK_CLIQUE, MAX_CLQ_SIZE, NB_ISETS_TEST, NB_INC_SUCCESS);
  
  fp_time = fopen("resultMaxCLQ", "a");
  fprintf(fp_time, 
	  "NewIncMaxCLQ %s %d %d %5.4f %d %d %d %d %d\n", 
	  argv[1], NB_NODE, NB_EDGE, 
	  ((float)NB_EDGE*2)/(NB_NODE*(NB_NODE-1)),
	  ((int)(endtime.ru_utime.tv_sec - starttime.ru_utime.tv_sec)), 
	  NB_BACK_CLIQUE, MAX_CLQ_SIZE, NB_ISETS_TEST, NB_INC_SUCCESS);
  fclose(fp_time);

  return TRUE;
}
