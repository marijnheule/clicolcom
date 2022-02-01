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

  int a, b;
  int size = 0;
  while (1) {
    int tmp = fscanf (graph, " e %i %i ", &a, &b);
    if (tmp == 0 || tmp == EOF) break;
    in[a-1] = 1; in[b-1] = 1;
    edges[size++] = a - 1;
    edges[size++] = b - 1; }

  int active = 0;
  for (int i = 0; i < nVertex; i++)
    if (in[i]) active++;

  int *verts = (int*) malloc (sizeof(int) * nVertex);

  for (int i = 0; i < nVertex; i++) verts[i] = i;

  size = 0;
  if (argc > 3) {
    FILE *order = fopen (argv[3], "r");
    int v;
    while (1) {
      int tmp = fscanf (order, " %i ", &v);
      if (tmp == 0 || tmp == EOF) break;
      if (v > 0)
        verts[size++] = v - 1; } }

  int extra = 0;
  if (argc > 3)
    extra = active * (colors - 1);

  printf ("p cnf %i %i\n", nVertex * colors, active + extra + nEdge * colors);

  for (int i = 0; i < nVertex; i++)
    if (in[i]) {
      for (int j = 1; j <= colors; j++)
        printf ("%i ", i * colors + j);
      printf ("0\n"); }

  if (argc > 3)
    for (int i = 0; i < nVertex; i++)
      if (in[i])
        for (int j = 2; j <= colors; j++) {
          for (int k = 0; k < i; k++)
            if (in[k]) printf ("%i ", verts[k] * colors + j - 1);
          printf ("-%i 0\n", verts[i] * colors + j); }

  for (int i = 0; i < nEdge; i++)
    for (int j = 1; j <= colors; j++)
      printf ("-%i -%i 0\n", edges[2*i] * colors + j, edges[2*i + 1] * colors + j);
}
