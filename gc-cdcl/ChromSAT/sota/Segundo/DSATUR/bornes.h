#include <float.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>

#include "graphe.h"
#include "cliquer/graph.h"
#include "cliquer/cliquer.h"
#include "cliquer/reorder.h"

//Calcul de bornes
int clique_cliquer_init(C_Graphe &GC, double cutoff, int *solution, int *size_solution);
