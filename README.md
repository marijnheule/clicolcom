# clicolcom
Code for the Cliques, Colors, and Communication project

Build using: ./build.sh

To compute the chromatic number of a graph run:

  ./solve.sh graph.edge

in which graph.edge is the graph in the DIMACS format.

The graphs in the DIMACS format start with a header,
which begins with "p edge " followed by the maximum
vertex index and the number of edges (seperated by a
space). Afterward, each line represents an edge and 
starts with "e " followed by the two vertices that
are connected by the edge (seperated by a space).
Note that all vertex labels should have a number
between 1 and the maximum vertex index. 

For example, a 5-cycle in DIMACS is:

p edge 5 5
e 1 2
e 2 3
e 3 4
e 4 5
e 5 1
