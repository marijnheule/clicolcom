#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <assert.h>

//#define DIRECT
#define SINZ

//#define WCNF

// argv[1] : graph in DIMACS format
// argv[2] : number of colors
// argv[3] : graph coloring (optional)
// argv[4] : offset (optional, default = 0)
// argv[5] : number of colorings to use

void atmostk (int *array, int size, int bound, int start) {
  // initilize the top right
  printf("-%i 0\n", start + (bound + 1) * size);

  // initialize the first row
  for (int j = 1; j <= size; j++) {
    printf("-%i %i 0\n", array[j-1], start + j); }

  // horizontal implications and increment implications
  for (int i = 0; i <= bound; i++)
    for (int j = 1; j < size; j++)
      printf("-%i %i 0\n", start + size * i + j, start + size * i + j + 1);

  // increment implications
  for (int i = 0; i < bound; i++)
    for (int j = 1; j < size; j++)
      printf("-%i -%i %i 0\n", array[j], start + size * i + j, start + size * (i+1) + j + 1);
}

void atleastk (int *array, int size, int bound, int start) {
  // make the lower left corner equal to the first representative
  printf("%i -%i 0\n", array[0], start + 1);
  if (bound > 1) printf("-%i 0\n", start + 1 + size);
  printf("%i 0\n", start + bound * size);

  // horizontal implications and increment implications
  for (int i = 0; i < bound; i++)
    for (int j = 1; j < size; j++) {
      printf("%i %i -%i 0\n", array[j], start + size * i + j, start + size * i + j + 1); }

  // diagonal implications
  for (int i = 0; i < bound - 1; i++)
    for (int j = 1; j < size; j++)
      printf("%i -%i 0\n", start + size * i + j, start + size * (i+1) + j + 1);
}

void subsetrec(int *in, int size, int max, int idx, int *out, int depth, int sign) {
  if (idx == max) {
    for (int j = 0; j < max; j++) {
      if (sign == -1) printf ("-");
      printf("%i ", out[j]); }
    printf("0\n");
    return; }

  if (depth >= size) return;

  out[idx] = in[depth];
  subsetrec (in, size, max, idx+1, out, depth+1, sign);
  subsetrec (in, size, max, idx  , out, depth+1, sign); }

int main (int argc, char** argv) {
  if (argc <= 2) {
    printf ("c run using ./maxclique GRAPH COLOR SOLUTION [offset]\n");
    exit (0); }

  FILE *graph  = fopen (argv[1], "r");
  int nColor = atoi (argv[2]);

  int nVertex, nEdge;
  fscanf (graph, " p edge %i %i ", &nVertex, &nEdge);

  int **adj;
  adj = (int **) malloc (sizeof(int *) * nVertex);
  for (int i = 0; i < nVertex; i++) {
    adj[i] = (int *) malloc (sizeof(int) * nVertex);
    for (int j = 0; j < nVertex; j++) adj[i][j] = 0; }

  int size = 0;
  while (1) {
    int a, b;
    int tmp = fscanf (graph, " e %i %i ", &a, &b);
    if (tmp == 0 || tmp == EOF) break;
    size++;
    adj[a-1][b-1] = 1;
    adj[b-1][a-1] = 1; }

  assert (size == nEdge);

  fclose (graph);

  int offset = 0;
  int numColorings = 1;
  if (argc > 4)
  {
    offset = atoi (argv[4]);
    numColorings = atoi (argv[5]);
  }

  int nVar = nVertex;

  if (argc <= 3) nVar += nVertex * nColor;

  if (argc > 3) {
#ifdef SINZ
    int numSinzVarsPerColoring = (nColor * (offset+2));
    if (offset) nVar += (numSinzVarsPerColoring * numColorings);
#endif
#ifdef DIRECT
    if (offset) nVar += nColor;
#endif
  }

  int nCls = nVertex * (nVertex - 1) / 2 - nEdge + nColor;
  if (argc <= 3) {
    nCls += (2*nColor - 1) * (nVertex - 1) + 3 - nColor;
    if (nColor == 1) nCls--; }

  if (argc > 3) {
#ifdef SINZ
    int numSinzClsPerColoring = (2*offset + 2) * (nColor - 1) + 2;
    if (offset > 0)
    {
      nCls += (numSinzClsPerColoring * numColorings);
    }
#endif
#ifdef DIRECT
    if (offset > 0) {
      int factor = nColor;
      for (int i = 1; i <= offset; i++)
        factor = factor * (nColor - i) / (i+1);
      nCls += factor; }
#endif
  }

#ifdef WCNF
  printf ("p wcnf %i %i %i\n", nVar, nCls, nVertex * (nVertex - 1) / 2 - nEdge);
#else
  printf ("p cnf %i %i\n", nVar, nCls);
#endif

  for (int i = 0; i < nVertex; i++)
    for (int j = i + 1; j < nVertex; j++)
      if (adj[i][j] == 0) {
#ifdef WCNF
        printf ("1 ");
#endif
        printf ("-%i -%i 0\n", i+1, j+1); }


  if (argc <= 3) {
    int *in  = (int *) malloc (sizeof (int) * nVertex);
    for (int i = 0; i < nVertex; i++) in [i] = i + 1;

    atleastk (in, nVertex, nColor, nVertex);
    exit (0);
  }

  if (argc > 3) {
    struct dirent *de;
    DIR *dr = opendir(argv[3]);
    if (dr == NULL)
    {
      return 0;
    }
    else
    {
      while ((de = readdir(dr)) != NULL)
      {
        FILE *sol = fopen (de->d_name, "r");
        int *in  = (int *) malloc (sizeof (int) * nColor);
        for (int i = 0; i < nColor; i++) in [i] = nVertex + i + 1 + (nColor * (numColorings - 1));

        int *color = (int*) malloc (sizeof(int) * nVertex);
        for (int i = 0; i < nVertex; i++) color[i] = 0;

        int lit;
        while (1) {
          int tmp = fscanf (sol, " %i ", &lit);
          if (tmp == 0 || tmp == EOF) break;
          if (lit > 0)
            color [(lit-1)/nColor] = ((lit - 1) % nColor) + 1; }

        fclose (sol);

        for (int c = 1; c <= nColor; c++) {
          if (offset) printf ("%i ", nVertex + c + (nColor * (numColorings - 1)));
#ifdef WCNF
          printf ("%i ", nVertex * (nVertex - 1) / 2 - nEdge);
#endif
          for (int i = 0; i < nVertex; i++)
            if (color[i] == c) printf ("%i ", i + 1);
          printf ("0\n"); }

        if (offset > 0) {
#ifdef SINZ
          atmostk (in, nColor, offset, nVertex + (nColor * (numColorings - 1)));
#endif
#ifdef DIRECT
          int *out = (int *) malloc (sizeof (int) * nColor);
          subsetrec (in, nColor, offset+1, 0, out, 0, -1);
#endif
        }
      }
      closedir(dr);
    }
  }
}
