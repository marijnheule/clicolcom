#include <stdio.h>
#include <stdlib.h>

// Given a graph G a parameter k, the code produces
// a formula that is satisfiable iff G is k-colorable.

// The graphs has an optional third input which is a
// vertex order. It will be used to construct a
// symmetry-breaking predicate

int main (int argc, char** argv) {
  FILE *graph;
  int colors, nVertex, nEdge;

  graph  = fopen (argv[1], "r");
  colors = atoi  (argv[2]);

  fscanf (graph, " p edge %i %i ", &nVertex, &nEdge);

  int *edges = (int*) malloc (sizeof(int) * 2*nEdge);
  int *in    = (int*) malloc (sizeof(int) * nVertex);
  for (int i = 0; i < nVertex; i++) in[i] = 0;

  int **adj;
  adj = (int **) malloc (sizeof(int *) * nVertex);
  for (int i = 0; i < nVertex; i++) {
    adj[i] = (int *) malloc (sizeof(int) * nVertex);
    for (int j = 0; j < nVertex; j++) adj[i][j] = 0; }

  int a, b;
  int size = 0;
  while (1) {
    int tmp = fscanf (graph, " e %i %i ", &a, &b);
    if (tmp == 0 || tmp == EOF) break;
    in[a-1] = 1; in[b-1] = 1;
    adj[a-1][b-1] = 1;
    adj[b-1][a-1] = 1;
    edges[size++] = a - 1;
    edges[size++] = b - 1; }

  int active = 0;
  for (int i = 0; i < nVertex; i++)
    if (in[i]) active++;

  int *verts = (int*) malloc (sizeof(int) * nVertex);

  for (int i = 0; i < nVertex; i++) verts[i] = i;

  size = 0;
  int clique = 0;
  int *fixed = (int*) malloc (sizeof(int) * nVertex);
  for (int i = 0; i < nVertex; i++) fixed[i] = 0;
  if (argc > 3) {
    FILE *order = fopen (argv[3], "r");
    int v;
    while (1) {
      int tmp = fscanf (order, " %i ", &v);
      if (tmp == 0 || tmp == EOF) break;
      if (v > 0)
        verts[size++] = v - 1; }

    for (int i = 0; i < nVertex; i++) {
      int flag = 1;
      for (int j = 0; j < i; j++)
        if (adj[verts[i]][verts[j]] == 0) flag = 0;
      if (flag) {
        clique++;
        fixed[verts[i]] = clique; }
      else break; } }
/*
  printf ("c clique size %i\n", clique);
  for (int i = 0; i < clique; i++)
    printf ("c %i %i\n", verts[i], fixed[verts[i]]);
*/
  int extra = 0;
  if (argc > 3 && size > clique)
    extra = (size - clique) * (colors - 1);
//    extra = (active - clique) * (colors - 1);

  int nCls = active + extra - clique;
  for (int i = 0; i < nEdge; i++) {
    int a = edges[2*i  ];
    int b = edges[2*i+1];
    if (!fixed[a] && !fixed[b]) nCls += colors; }

  printf ("p cnf %i %i\n", nVertex * colors, nCls);

  int* domain = (int*) malloc (sizeof (int) * (colors + 1));
  for (int i = 0; i < nVertex; i++)
    if (in[i]) {
      if (!fixed[i]) { // TODO: also fix clique?
        for (int j = 1; j <= colors; j++)
          domain[j] = 1;
        for (int j = 0; j < nVertex; j++)
          if (adj[j][i]) domain[fixed[j]] = 0;
        for (int j = 1; j <= colors; j++)
          if (domain[j]) printf ("%i ", i * colors + j);
        printf ("0\n"); } }

  if (argc > 3 && size > clique)
    for (int i = clique; i < size; i++)
      if (in[i])
        for (int j = 2; j <= colors; j++) {
          for (int k = clique; k < i; k++)
            if (in[k]) printf ("%i ", verts[k] * colors + j - 1);
          printf ("-%i 0\n", verts[i] * colors + j); }

  for (int i = 0; i < nEdge; i++) {
    int a = edges[2*i  ];
    int b = edges[2*i+1];
    if (!fixed[a] && !fixed[b])
      for (int j = 1; j <= colors; j++)
        printf ("-%i -%i 0\n", a * colors + j, b * colors + j); }
}
