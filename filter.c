#include <stdio.h>
#include <stdlib.h>

// Given a graph and a set of vertices, the
// code produces the induces subgraph

int main (int argc, char** argv) {
  FILE *graph, *core;
  int nVertex, nEdge;

  graph = fopen (argv[1], "r");
  core  = fopen (argv[2], "r");

  fscanf (graph, " p edge %i %i ", &nVertex, &nEdge);

  int *in = (int*) malloc (sizeof(int) * nVertex);

  int v;
  while (1) {
    int tmp = fscanf (core, " %i ", &v);
    if (tmp == 0 || tmp == EOF) break;
    in[v-1] = 1; }

  int a, b;
  int size = 0;
  while (1) {
    int tmp = fscanf (graph, " e %i %i ", &a, &b);
    if (tmp == 0 || tmp == EOF) break;
    if (in[a-1] && in[b-1])
      printf ("e %i %i\n", a, b); }
}
