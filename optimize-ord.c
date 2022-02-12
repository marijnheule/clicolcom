#include <stdio.h>
#include <stdlib.h>

// argv[1] : graph in DIMACS format
// argv[2] : partial ordering as a sequence of vertex numbers

int main (int argc, char** argv) {
  FILE *graph;
  int nVertex, nEdge;
  int details = 0;

  graph  = fopen (argv[1], "r");

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

  fclose (graph);

  int *ord = (int*) malloc (sizeof(int) * nVertex);
  for (int i = 0; i < nVertex; i++) ord[i] = -1;

  FILE *order = fopen (argv[2], "r");
  int v, c = 0;
  while (1) {
    int tmp = fscanf (order, " %i ", &v);
    if (tmp == 0 || tmp == EOF) break;
    if (v < nVertex) {
      ord[c++] = v - 1;
      in[v-1] = 0; } }

  fclose (order);

  if (argc > 3) {
    if (argv[3][0] == '-' && argv[3][1] == 'v') details = 1; }

  for (int i = 0; i < nVertex; i++) {
    if (in[i]) {
      ord[c++] = i;
      in[i] = 0; } }

  int *mask = (int*) malloc (sizeof(int) * (nVertex));
  for (int i = 0; i <= nVertex; i++) mask[i] = 0;

  int clique = 0;
  for (int i = 0; i < nVertex; i++) {
    int count = 0;
    for (int j = 0; j < size; j += 2) {
      if (edges[j  ] == ord[i]) { if (mask[edges[j+1]]) count++; }
      if (edges[j+1] == ord[i]) { if (mask[edges[j  ]]) count++; }
    }

    if (i <= 3 * clique)
    if (count != i) {
      int swap = i;
      int best = count;
      for (int k = i; k < nVertex; k++) {
        int count = 0;
        for (int j = 0; j < size; j += 2) {
          if (edges[j  ] == ord[k]) { if (mask[edges[j+1]]) count++; }
          if (edges[j+1] == ord[k]) { if (mask[edges[j  ]]) count++; } }
        if (count > best) { swap = k; best = count; } }
      if (swap != i) {
        int o = ord[i];
        ord[i] = ord[swap];
        ord[swap] = o;
      }
    }

    mask[ord[i]] = 1;
    count = 0;
    for (int j = 0; j < size; j += 2) {
      if (edges[j  ] == ord[i]) { if (mask[edges[j+1]]) count++; }
      if (edges[j+1] == ord[i]) { if (mask[edges[j  ]]) count++; }
    }
    if (count == i) clique++;
    if (details)
      printf ("c vertex %i is connected to %i predecessors\n", ord[i] + 1, count);

    printf ("%i\n", ord[i] + 1);

  }
  if (details)
    printf ("c top clique of size %i\n", clique);

}
