#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CONVERTER_H_
#define CONVERTER_H_

typedef struct edge Edge;

struct edge{
	int a,b;
};

FILE * daOutput;
FILE * daInput;

/* The following code is copied from micheal trick's implemantation of the DSATUR algorithm */
/* mat.gsia.cmu.edu/COLOR/solvers/trick.c */

#define MAX_NODE 6000
#define TRUE 1
#define FALSE 0

int adj[MAX_NODE][MAX_NODE];
int BestColoring;
int ColorClass[MAX_NODE];
int prob_count;
int Order[MAX_NODE];
int Handled[MAX_NODE];
int ColorAdj[MAX_NODE][MAX_NODE];
int ColorCount[MAX_NODE];
int valid[MAX_NODE],clique[MAX_NODE];
int neighbours[MAX_NODE][MAX_NODE];
int lb;
int best_clique;
int num_prob,max_prob;

int p_mergeliterals;
int p_clique;
int nr_of_nodes;
int nr_of_edges;
int nr_of_colours;

int nr_of_discarded_nodes;

Edge * edges;
int * degree_info;

void init_vars(int nrnodes, int nredges, int nrcolours, FILE * edgeFile);

void find_clique();

void free_vars();

void print_header();

void print_at_least_one();

void print_at_most_one();

void print_conflict();

void print_e();

void print_clique();

void print_dummies();

#endif /*CONVERTER_H_*/
