/******************************************************************************/
/*                Graph coloring using techologies of Satz                    */
/*                Author: Chu-Min LI                                          */
/*                Copyright MIS, University of Picardie Jules Verne, France   */
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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define NONE2 -2

#define WEIGHT 4
#define WEIGHT1 25
#define WEIGHT2 5
#define WEIGHT3 1
#define T 10

/* the tables of nodes and edges are statically allocated. Modify the
   parameters tab_node_size and tab_edge_size before compilation if
   necessary */
#define tab_node_size 11000
#define tab_edge_size 1100000
#define max_nb_value 500
#define pop(stack) stack[--stack##_fill_pointer]
#define push(item, stack) stack[stack##_fill_pointer++] = item

#define PASSIVE 0
#define ACTIVE 1

enum algo_type { TOP_DOWN, BOTTOM_UP, BINARY };

int *node_neibors[tab_node_size];
int node_value[tab_node_size];
int node_state[tab_node_size];
int *node_value_state[tab_node_size];
int node_nb_value[tab_node_size];
int node_nb_symmetric_value[tab_node_size];
int node_impact[tab_node_size];
int nb_neibors[tab_node_size];
int static_nb_neibors[tab_node_size];
int branching_node[tab_node_size];
int NODE_STACK[tab_node_size];
int NODE_STACK_fill_pointer = 0;
int CHANGE_STACK[2 * tab_node_size * max_nb_value];
int CHANGE_STACK_fill_pointer = 0;
int UNITNODE_STACK[tab_node_size];
int UNITNODE_STACK_fill_pointer = 0;

int CANDIDATE_NODE_STACK[tab_node_size];
int CANDIDATE_NODE_STACK_fill_pointer = 0;

int saved_change_stack[tab_node_size];

int edge[tab_edge_size][2];
int symmetry[max_nb_value];
int broken_symmetry_stack[max_nb_value];
int broken_symmetry_stack_fill_pointer = 0;
int saved_broken_symmetry_stack[tab_node_size];

int NB_VALUE, NB_NODE, NB_EDGE, H;

int NB_UNIT = 0, NB_BACK = 0, NB_BRANCHE = 0;

int lower_bound;
int upper_bound;
int policy;

// std::mt19937 random_generator;

int edge_redundant(int node1, int node2, int nb_edges) { /* node1<node2 */
  int i;

  for (i = 0; i < nb_edges; i++) {
    if (node1 == edge[i][0] && node2 == edge[i][1])
      return TRUE;
  }
  return FALSE;
}

int FORMAT = 1;

void free_graph_instance() {
  int i;
  for (i = 1; i <= NB_NODE; i++) {
    free(node_neibors[i]);
    free(node_value_state[i]);
  }
}

my_type build_simple_graph_instance(char *input_file) {
  printf("read file %s\n", input_file);

  FILE *fp_in = fopen(input_file, "r");
  char ch, word2[WORD_LENGTH];
  int i, j, e, node1, node2;
  if (fp_in == NULL)
    return FALSE;
  if (FORMAT == 1) {
    fscanf(fp_in, "%c", &ch);
    while (ch != 'p') {
      while (ch != '\n')
        fscanf(fp_in, "%c", &ch);
      fscanf(fp_in, "%c", &ch);
    }
    fscanf(fp_in, "%s%d%d", word2, &NB_NODE, &NB_EDGE);
  } else
    fscanf(fp_in, "%d%d", &NB_NODE, &NB_EDGE);

  // printf(" #nodes=%i, #edges=%i\n", NB_NODE, NB_EDGE);

  for (i = 1; i <= NB_NODE; i++) {
    nb_neibors[i] = 0;
  }

  for (i = 0; i < NB_EDGE; i++) {
    if (FORMAT == 1) {
      fscanf(fp_in, "%s%d%d", word2, &edge[i][0], &edge[i][1]);
    } else
      fscanf(fp_in, "%d%d", &edge[i][0], &edge[i][1]);

    // if((i%1000) == 0)
    // 	printf(" %i-th edge=%i,%i\n", i, edge[i][0], edge[i][1]);

    if (edge[i][0] == edge[i][1]) {
      i--;
      NB_EDGE--;
      // printf("auto edge %d over %d\n", i--, NB_EDGE--);
    } else {
      if (edge[i][0] > edge[i][1]) {
        e = edge[i][1];
        edge[i][1] = edge[i][0];
        edge[i][0] = e;
      }
      if (edge_redundant(edge[i][0], edge[i][1], i) == TRUE) {
        i--;
        NB_EDGE--;
        // printf("edge redundant %d over %d", i--, NB_EDGE--);
      } else {
        nb_neibors[edge[i][0]]++;
        nb_neibors[edge[i][1]]++;
      }
    }
  }
  fclose(fp_in);

  // printf("here\n");

  for (i = 1; i <= NB_NODE; i++) {
    static_nb_neibors[i] = nb_neibors[i];
    node_neibors[i] = (int *)malloc((nb_neibors[i] + 1) * sizeof(int));
    node_neibors[i][nb_neibors[i]] = NONE;
    nb_neibors[i] = 0;
    node_nb_value[i] = NB_VALUE;
    node_nb_symmetric_value[i] = NB_VALUE;
    node_state[i] = ACTIVE;
    node_value_state[i] = (int *)malloc(NB_VALUE * sizeof(int));
    for (j = 0; j < NB_VALUE; j++)
      node_value_state[i][j] = ACTIVE;
  }
  for (i = 0; i < NB_EDGE; i++) {
    node1 = edge[i][0];
    node2 = edge[i][1];
    node_neibors[node1][nb_neibors[node1]++] = node2;
    node_neibors[node2][nb_neibors[node2]++] = node1;
  }
  // printf("there\n");

  return TRUE;
}

int assign_value(int node, int value) {
  int i, *neibors, neibor;

  neibors = node_neibors[node];
  for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
    if (node_state[neibor] == ACTIVE &&
        node_value_state[neibor][value] == ACTIVE) {
      node_value_state[neibor][value] = PASSIVE;
      push(neibor, CHANGE_STACK);
      push(value, CHANGE_STACK);
      node_nb_value[neibor]--;
      if (node_nb_value[neibor] == 0)
        return NONE;
      else if (node_nb_value[neibor] == 1)
        push(neibor, UNITNODE_STACK);
    }
  }
  node_value[node] = value;
  node_state[node] = PASSIVE;
  push(node, NODE_STACK);

  if (symmetry[value] == TRUE) {
    push(value, broken_symmetry_stack);
    symmetry[value] = FALSE;
    /*
    for(i=1; i<=NB_NODE; i++) {
      if (node_state[i]==ACTIVE && node_value_state[i][value]==ACTIVE)
        node_nb_symmetric_value[i]--;
    }
    */
  }
  return TRUE;
}

int get_active_value(int node) {
  int i;
  for (i = 0; i < NB_VALUE; i++) {
    if (node_value_state[node][i] == ACTIVE) {
      return i;
    }
  }
  return NONE;
}

int unitnode_propagation() {
  int unitnode_position, unitnode, value;
  for (unitnode_position = 0; unitnode_position < UNITNODE_STACK_fill_pointer;
       unitnode_position++) {
    unitnode = UNITNODE_STACK[unitnode_position];
    if (node_state[unitnode] == ACTIVE) {
      NB_UNIT++;
      value = get_active_value(unitnode);
      if (value == NONE)
        printf("error...\n");
      if (assign_value(unitnode, value) == NONE)
        return NONE;
      branching_node[unitnode] = FALSE;
    }
  }
  UNITNODE_STACK_fill_pointer = 0;
  return TRUE;
}

void remove_symmetry(int node, int value) {
  int i;
  for (i = 0; i < NB_VALUE; i++) {
    if (i != value && (node_value_state[node][i] == ACTIVE) &&
        symmetry[i] == TRUE) {
      node_value_state[node][i] = PASSIVE;
      push(node, CHANGE_STACK);
      push(i, CHANGE_STACK);
      node_nb_value[node]--;
    }
  }
}

int backtracking() {
  int node, node1, value, i;

  NB_BACK++;
  UNITNODE_STACK_fill_pointer = 0;
  do {
    node = pop(NODE_STACK);
    if (branching_node[node] == FALSE)
      node_state[node] = ACTIVE;
    else if (node_nb_value[node] == 1)
      node_state[node] = ACTIVE;
    else {
      for (i = saved_change_stack[node]; i < CHANGE_STACK_fill_pointer;
           i = i + 2) {
        node1 = CHANGE_STACK[i];
        node_value_state[node1][CHANGE_STACK[i + 1]] = ACTIVE;
        node_nb_value[node1]++;
      }
      CHANGE_STACK_fill_pointer = saved_change_stack[node];
      for (i = saved_broken_symmetry_stack[node];
           i < broken_symmetry_stack_fill_pointer; i++) {
        value = broken_symmetry_stack[i];
        symmetry[value] = TRUE;
        /*
        for(node1=1; node1<=NB_NODE; node1++) {
          if (node_state[node1]==ACTIVE &&
        node_value_state[node1][value]==ACTIVE)
            node_nb_symmetric_value[node1]++;
        }
        */
      }
      broken_symmetry_stack_fill_pointer = saved_broken_symmetry_stack[node];
      node_nb_value[node]--;
      node_value_state[node][node_value[node]] = PASSIVE;
      push(node, CHANGE_STACK);
      push(node_value[node], CHANGE_STACK);
      value = get_active_value(node);
      if (symmetry[value] == TRUE)
        remove_symmetry(node, value);
      saved_change_stack[node] = CHANGE_STACK_fill_pointer;
      if ((value == NONE) || assign_value(node, value) == NONE) {
        printf("erreur....");
        exit(1);
        return NONE;
      }
      return TRUE;
    }
  } while (NODE_STACK_fill_pointer > 0);
  return FALSE;
}

int verify_solution() {
  int i;

  for (i = 0; i < NB_EDGE; i++) {
    if (node_value[edge[i][0]] == node_value[edge[i][1]])
      return FALSE;
  }
  return TRUE;
}

int node_nb_eliminated_value[tab_node_size];
void check_consistency() {
  int i, j, nb_passive_node = 0, nb_passive_value;

  // all passive nodes should be in NODE_STACK and all nodes
  // in NODE_STACK should be passive
  for (i = 1; i <= NB_NODE; i++) {
    if (node_state[i] == PASSIVE)
      nb_passive_node++;
  }
  for (i = 0; i < NODE_STACK_fill_pointer; i++)
    if (node_state[NODE_STACK[i]] == ACTIVE)
      printf("erreur1...\n");
  if (nb_passive_node != NODE_STACK_fill_pointer)
    printf("erreur1...\n");

  // all eliminated value should be passive and should be in CHANGE_STACK
  // the number of eliminated values plus the number of remaining values
  // of a node should be the total number of values (NB_VALUE)
  for (i = 1; i <= NB_NODE; i++) {
    node_nb_eliminated_value[i] = 0;
  }
  for (i = 0; i < CHANGE_STACK_fill_pointer; i = i + 2) {
    if (node_value_state[CHANGE_STACK[i]][CHANGE_STACK[i + 1]] == ACTIVE)
      printf("erreur2...\n");
    node_nb_eliminated_value[CHANGE_STACK[i]]++;
  }
  for (i = 1; i <= NB_NODE; i++) {
    nb_passive_value = 0;
    for (j = 0; j < NB_VALUE; j++) {
      if (node_value_state[i][j] == PASSIVE)
        nb_passive_value++;
    }
    if ((nb_passive_value != node_nb_eliminated_value[i]) ||
        nb_passive_value + node_nb_value[i] != NB_VALUE)
      printf("erreur2...\n");
  }
}

int choose_node_dom() {
  int node, i, min_nb_value;
  min_nb_value = NB_VALUE + 1;
  node = NONE;
  NB_BRANCHE++;
  for (i = 1; i <= NB_NODE; i++) {
    if (node_state[i] == ACTIVE) {
      if (node_nb_value[i] < min_nb_value) {
        min_nb_value = node_nb_value[i];
        node = i;
      }
    }
  }
  if (node == NONE) {
    printf("erreur...");
    exit(1);
  }
  return node;
}

void choose_candidates() {
  int i, min_nb_value;
  min_nb_value = NB_VALUE + 1;
  CANDIDATE_NODE_STACK_fill_pointer = 0;
  for (i = 1; i <= NB_NODE; i++) {
    if (node_state[i] == ACTIVE) {
      if (node_nb_value[i] < min_nb_value) {
        min_nb_value = node_nb_value[i];
        CANDIDATE_NODE_STACK_fill_pointer = 0;
        push(i, CANDIDATE_NODE_STACK);
      } else if (node_nb_value[i] == min_nb_value)
        push(i, CANDIDATE_NODE_STACK);
    }
  }
}

int choose_by_degree() {
  int max_degree, nb_degree, *neibors, neibor, chosen_node = FALSE, node, i, nb_best;
  max_degree = -1;
	nb_best = 0;
  for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
    node = CANDIDATE_NODE_STACK[i];
    nb_degree = 0;
    neibors = node_neibors[node];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (node_state[neibor] == ACTIVE)
        nb_degree++;
    }
    if (nb_degree > max_degree) {
      max_degree = nb_degree;
			chosen_node = node;
			nb_best = 1;
    } else if(nb_degree == max_degree) {
    	++nb_best;
			if(rand() % nb_best == 0) {
					chosen_node = node;
			}
    }
  }
  return chosen_node;
}

int choose_node_dom_deg() {
  int node, i, min_nb_value, chosen_node;
  NB_BRANCHE++;
  choose_candidates();
  if (CANDIDATE_NODE_STACK_fill_pointer == 0) {
    printf("erreur...");
    exit(1);
  } else if (CANDIDATE_NODE_STACK_fill_pointer == 1)
    chosen_node = CANDIDATE_NODE_STACK[0];
  else
    chosen_node = choose_by_degree();
  return chosen_node;
}

int choose_by_neibor_degree() {
  int max_degree, nb_degree, *neibors, neibor, chosen_node = FALSE, node, i,
                                               min_nb_value;
  max_degree = -1;
  min_nb_value = node_nb_value[CANDIDATE_NODE_STACK[0]];
  for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
    node = CANDIDATE_NODE_STACK[i];
    nb_degree = 0;
    neibors = node_neibors[node];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (node_state[neibor] == ACTIVE)
        switch (node_nb_value[neibor] - min_nb_value) {
        case 0:
          nb_degree = nb_degree + 3;
          break;
        case 1:
          nb_degree = nb_degree + 2;
          break;
        case 2:
          nb_degree = nb_degree + 1;
          break;
        default:
          nb_degree++;
        }
    }
    if (nb_degree > max_degree) {
      max_degree = nb_degree;
      chosen_node = node;
    }
  }
  return chosen_node;
}

int choose_by_neibor_degree2() {
  int max_degree, nb_degree, *neibors, neibor,
      chosen_node = FALSE, node, i, min_nb_value, start, my_neibor, *my_neibors;
  max_degree = -1;
  start = 0;
  min_nb_value = node_nb_value[CANDIDATE_NODE_STACK[0]];
  for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
    node = CANDIDATE_NODE_STACK[i];
    nb_degree = 0;
    neibors = node_neibors[node];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (node_state[neibor] == ACTIVE)
        switch (node_nb_value[neibor] - min_nb_value) {
        case 0:
          nb_degree = nb_degree + 3;
          break;
        case 1:
          nb_degree = nb_degree + 2;
          break;
        case 2:
          nb_degree = nb_degree + 1;
          break;
        default:
          nb_degree++;
        }
    }
    if (nb_degree > max_degree) {
      max_degree = nb_degree;
      CANDIDATE_NODE_STACK[0] = node;
      start = 1;
    } else if (nb_degree == max_degree) {
      CANDIDATE_NODE_STACK[start++] = node;
    }
  }
  if (start == 1)
    return CANDIDATE_NODE_STACK[0];
  max_degree = -1;
  for (i = 0; i < start; i++) {
    node = CANDIDATE_NODE_STACK[i];
    nb_degree = 0;
    neibors = node_neibors[node];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (node_state[neibor] == ACTIVE) {
        my_neibors = node_neibors[neibor];
        for (my_neibor = *my_neibors; my_neibor != NONE;
             my_neibor = *(++my_neibors)) {
          if (node_state[my_neibor] == ACTIVE)
            nb_degree++;
        }
      }
    }
    if (nb_degree > max_degree) {
      max_degree = nb_degree;
      chosen_node = node;
    }
  }
  return chosen_node;
}

int choose_by_neibor_degree1() {
  int max_degree, nb_degree, *neibors, neibor,
      chosen_node = FALSE, node, i, min_nb_value, nb_min_neibor,
      neibor_min_nb_value, min_neibor_min_nb_value, max_nb_min_neibor;
  max_degree = -1;
  nb_min_neibor = NB_NODE;
  min_nb_value = node_nb_value[CANDIDATE_NODE_STACK[0]];
  for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
    node = CANDIDATE_NODE_STACK[i];
    nb_degree = 0; // total_nb_value=0;
    neibors = node_neibors[node];
    neibor_min_nb_value = NB_VALUE + 1;
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (node_state[neibor] == ACTIVE) {
        nb_degree++; // total_nb_value=total_nb_value+node_nb_value[neibor];
        if (node_nb_value[neibor] < neibor_min_nb_value) {
          neibor_min_nb_value = node_nb_value[neibor];
          nb_min_neibor = 1;
        } else if (node_nb_value[neibor] == neibor_min_nb_value)
          nb_min_neibor++;
      }
    }
    if (nb_degree > max_degree) {
      max_degree = nb_degree;
      chosen_node = node;
      min_neibor_min_nb_value = neibor_min_nb_value;
      max_nb_min_neibor = nb_min_neibor;
    } else if ((nb_degree == max_degree) &&
               (min_neibor_min_nb_value > neibor_min_nb_value)) {
      chosen_node = node;
      min_neibor_min_nb_value = neibor_min_nb_value;
      max_nb_min_neibor = nb_min_neibor;
    }
    /*
    else if ((nb_degree==max_degree) && (
    min_neibor_min_nb_value==neibor_min_nb_value) &&
             (max_nb_min_neibor<nb_min_neibor)) {
      chosen_node=node;
      max_nb_min_neibor=nb_min_neibor;
    }
    */
  }
  return chosen_node;
}

int my_assign_value(int node, int value) {
  int i, *neibors, neibor;

  neibors = node_neibors[node];
  for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
    if (node_state[neibor] == ACTIVE &&
        node_value_state[neibor][value] == ACTIVE) {
      node_value_state[neibor][value] = PASSIVE;
      push(neibor, CHANGE_STACK);
      push(value, CHANGE_STACK);
      node_nb_value[neibor]--;
      if (node_nb_value[neibor] == 0)
        return NONE;
      else if (node_nb_value[neibor] == 1)
        push(neibor, UNITNODE_STACK);
    }
  }
  node_value[node] = value;
  node_state[node] = PASSIVE;
  push(node, NODE_STACK);
  return TRUE;
}

int my_unitnode_propagation() {
  int unitnode_position, unitnode, value;
  for (unitnode_position = 0; unitnode_position < UNITNODE_STACK_fill_pointer;
       unitnode_position++) {
    unitnode = UNITNODE_STACK[unitnode_position];
    if (node_state[unitnode] == ACTIVE) {
      NB_UNIT++;
      value = get_active_value(unitnode);
      if (value == NONE)
        printf("error...\n");
      if (my_assign_value(unitnode, value) == NONE)
        return NONE;
    }
  }
  UNITNODE_STACK_fill_pointer = 0;
  return TRUE;
}

int my_assign_value_and_propagate(int node, int value) {
  if (my_assign_value(node, value) == NONE) {
    printf("erreur\n");
    exit(1);
  }
  if (my_unitnode_propagation() == NONE)
    return NONE;
  return CHANGE_STACK_fill_pointer;
}

void reset_context(int saved_node_stack_fill_pointer,
                   int saved_change_stack_fill_pointer) {
  int i, node;
  for (i = saved_node_stack_fill_pointer; i < NODE_STACK_fill_pointer; i++)
    node_state[NODE_STACK[i]] = ACTIVE;
  NODE_STACK_fill_pointer = saved_node_stack_fill_pointer;
  for (i = saved_change_stack_fill_pointer; i < CHANGE_STACK_fill_pointer;
       i = i + 2) {
    node = CHANGE_STACK[i];
    node_value_state[node][CHANGE_STACK[i + 1]] = ACTIVE;
    node_nb_value[node]++;
  }
  CHANGE_STACK_fill_pointer = saved_change_stack_fill_pointer;
  UNITNODE_STACK_fill_pointer = 0;
}

int examine_candidate_value(int node, int value) {
  int impact, saved_node_stack_fill_pointer, saved_change_stack_fill_pointer,
      saved_broken_symmetry_stack_fill_pointer;
  saved_node_stack_fill_pointer = NODE_STACK_fill_pointer;
  saved_change_stack_fill_pointer = CHANGE_STACK_fill_pointer;
  saved_broken_symmetry_stack_fill_pointer = broken_symmetry_stack_fill_pointer;
  impact = my_assign_value_and_propagate(node, value);
  reset_context(saved_node_stack_fill_pointer, saved_change_stack_fill_pointer);
  if (impact == NONE) {
    node_nb_value[node]--;
    node_value_state[node][value] = PASSIVE;
    push(node, CHANGE_STACK);
    push(value, CHANGE_STACK);
    if (node_nb_value[node] == 0)
      return NONE;
    else if (node_nb_value[node] == 1) {
      push(node, UNITNODE_STACK);
      if (unitnode_propagation() == NONE)
        return NONE;
    }
    return FALSE;
  } else if (impact == saved_change_stack_fill_pointer) {
    assign_value(node, value);
    branching_node[node] = FALSE;
  } else if (impact > saved_change_stack_fill_pointer)
    node_impact[node] =
        node_impact[node] * (impact - saved_change_stack_fill_pointer);
  return impact - saved_change_stack_fill_pointer;
}

int examine_candidates() {
  int i, j, node, impact, symmetry_impact = NONE2;
  for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
    node = CANDIDATE_NODE_STACK[i];
    if (node_state[node] == ACTIVE) {
      node_impact[node] = 1;
      for (j = 0; j < NB_VALUE; j++) {
        if (node_value_state[node][j] == ACTIVE) {
          if (symmetry[j] == FALSE || symmetry_impact == NONE2) {
            impact = examine_candidate_value(node, j);
            if (impact == NONE)
              return NONE;
            else if (node_state[node] == PASSIVE)
              break;
            if (symmetry[j] == TRUE)
              symmetry_impact = impact;
          } else {
            if (symmetry_impact ==
                FALSE) { // symmetric value has led to a contradiction
              node_nb_value[node]--;
              node_value_state[node][j] = PASSIVE;
              push(node, CHANGE_STACK);
              push(j, CHANGE_STACK);
              if (node_nb_value[node] == 0)
                return NONE;
              else if (node_nb_value[node] == 1) {
                push(node, UNITNODE_STACK);
                if (unitnode_propagation() == NONE)
                  return NONE;
              }
            } else
              node_impact[node] = node_impact[node] * symmetry_impact;
          }
        }
      }
    }
  }
  return TRUE;
}

int choose_node_up() {
  int i, node, min_nb_value, max_impact, chosen_node = FALSE;
  NB_BRANCHE++;
  do {
    choose_candidates();
    /*  if (CANDIDATE_NODE_STACK_fill_pointer==0)
      printf("erreuraa...");
      else */
    if (CANDIDATE_NODE_STACK_fill_pointer == 1)
      chosen_node = CANDIDATE_NODE_STACK[0];
    else {
      min_nb_value = node_nb_value[CANDIDATE_NODE_STACK[0]];
      if (min_nb_value == 1) {
        printf("erreur\n");
        exit(1);
      } else if (min_nb_value == 2) {
        if (examine_candidates() == NONE)
          return NONE;
        max_impact = -1;
        for (i = 0; i < CANDIDATE_NODE_STACK_fill_pointer; i++) {
          node = CANDIDATE_NODE_STACK[i];
          if (node_state[node] == ACTIVE && node_impact[node] > max_impact) {
            max_impact = node_impact[node];
            chosen_node = node;
          }
        }
        // chosen_node=choose_by_degree();
      } else
        // chosen_node=choose_by_neibor_degree();
        chosen_node = choose_by_degree();
    }
  } while (chosen_node == FALSE && CANDIDATE_NODE_STACK_fill_pointer > 0);
  return chosen_node;
}

int search_by_up() {
  int node, value, i;
  if (unitnode_propagation() == NONE)
    return FALSE;
  do {
    //   check_consistency();
    if (NODE_STACK_fill_pointer < NB_NODE) {
      node = choose_node_up();
      if (node == NONE)
        backtracking();
      else if (node == FALSE) {
        if (NODE_STACK_fill_pointer == NB_NODE)
          return TRUE;
        return FALSE;
      } else {
        value = get_active_value(node);
        if (symmetry[value] == TRUE)
          remove_symmetry(node, value);
        branching_node[node] = TRUE;
        saved_change_stack[node] = CHANGE_STACK_fill_pointer;
        saved_broken_symmetry_stack[node] = broken_symmetry_stack_fill_pointer;
        if ((value == NONE) || assign_value(node, value) == NONE) {
          printf("erreur....");
          exit(1);
          return NONE;
        }
      }
      while ((NODE_STACK_fill_pointer > 0) && (unitnode_propagation() == NONE))
        backtracking();
    } else
      return TRUE;
  } while (NODE_STACK_fill_pointer != 0);
  return FALSE;
}

int search_by_dom_deg() {
  int node, value;
  if (unitnode_propagation() == NONE)
    return FALSE;
  do {
    // check_consistency();
    if (NODE_STACK_fill_pointer < NB_NODE) {
      node = choose_node_dom_deg();
      value = get_active_value(node);
      if (symmetry[value] == TRUE)
        remove_symmetry(node, value);
      branching_node[node] = TRUE;
      saved_change_stack[node] = CHANGE_STACK_fill_pointer;
      saved_broken_symmetry_stack[node] = broken_symmetry_stack_fill_pointer;
      if ((value == NONE) || assign_value(node, value) == NONE) {
        printf("erreur....");
        exit(1);
        return NONE;
      }
      while ((NODE_STACK_fill_pointer > 0) && (unitnode_propagation() == NONE))
        backtracking();
    } else
      return TRUE;
  } while (NODE_STACK_fill_pointer != 0);
  return FALSE;
}

int search_by_dom() {
  int node, value;
  if (unitnode_propagation() == NONE)
    return FALSE;
  do {
    // check_consistency();
    if (NODE_STACK_fill_pointer < NB_NODE) {
      node = choose_node_dom();
      value = get_active_value(node);
      if (symmetry[value] == TRUE)
        remove_symmetry(node, value);
      branching_node[node] = TRUE;
      saved_change_stack[node] = CHANGE_STACK_fill_pointer;
      saved_broken_symmetry_stack[node] = broken_symmetry_stack_fill_pointer;
      if ((value == NONE) || assign_value(node, value) == NONE) {
        printf("erreur....");
        exit(1);
        return NONE;
      }
      while ((NODE_STACK_fill_pointer > 0) && (unitnode_propagation() == NONE))
        backtracking();
    } else
      return TRUE;
  } while (NODE_STACK_fill_pointer != 0);
  return FALSE;
}

void init() {
  int i;
  NB_BACK = 0;
  NB_UNIT = 0;
  NB_BRANCHE = 0;
  NODE_STACK_fill_pointer = 0;
  CHANGE_STACK_fill_pointer = 0;
  UNITNODE_STACK_fill_pointer = 0;
  for (i = 0; i < NB_VALUE; i++)
    symmetry[i] = TRUE;
  broken_symmetry_stack_fill_pointer = 0;
}

void print_solution(char *instance, int NB_VALUE) {
  char sol[WORD_LENGTH];
  FILE *fp_sol;
  int i;
  sprintf(sol, "solution_for_%s_%d_color", instance, NB_VALUE);
  printf("solution verified is saved in %s file\n", sol);
  fp_sol = fopen(sol, "w");
  for (i = 1; i <= NB_NODE; i++) {
    printf("%d ", node_value[i] + 1);
    fprintf(fp_sol, "%d ", node_value[i] + 1);
  }
  printf("\n");
  fprintf(fp_sol, "\n");
  fclose(fp_sol);
}

void scanone(int argc, char *argv[], int i, int *varptr) {
  if (i >= argc || sscanf(argv[i], "%i", varptr) != 1) {
    fprintf(stderr, "Bad argument %s\n", i < argc ? argv[i] : argv[argc - 1]);
    exit(-1);
  }
}

int HELP_FLAG = FALSE;
char *INPUT_FILE;

void parse_parameters(int argc, char *argv[]) {
  int i, temp, j;

  lower_bound = 0;
  upper_bound = max_nb_value;
  policy = TOP_DOWN;
	int seed = 12345;

  if (argc < 2)
    HELP_FLAG = TRUE;
  else
    for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-nbColors") == 0)
        scanone(argc, argv, ++i, &NB_VALUE);
      else if (strcmp(argv[i], "-f") == 0)
        scanone(argc, argv, ++i, &FORMAT);
      else if (strcmp(argv[i], "-l") == 0)
        scanone(argc, argv, ++i, &lower_bound);
      else if (strcmp(argv[i], "-u") == 0)
        scanone(argc, argv, ++i, &upper_bound);
      else if (strcmp(argv[i], "-p") == 0)
        scanone(argc, argv, ++i, &policy);
      else if (strcmp(argv[i], "-s") == 0)
        scanone(argc, argv, ++i, &seed);
      else if (strcmp(argv[i], "-help") == 0)
        HELP_FLAG = TRUE;
      else
        INPUT_FILE = argv[i];
    }
		
		printf("random seed = %d\n", seed);
		srand(seed);
		
}

char *filename(char *input) {
  char c, *input1;
  int nb, nb1;
  input1 = input;
  nb = 0;
  for (c = *input1; c != '\0'; c = *(input1 += 1))
    if (c == '/')
      nb++;
  input1 = input;
  nb1 = 0;
  for (c = *input1; c != '\0'; c = *(input1 += 1)) {
    if (c == '/')
      nb1++;
    if (nb == nb1)
      return input1 + 1;
  }
  return input1;
}

struct tms *a_tms;

int solve(long begintime) {
  int i, result;
  long endtime, mess;
  FILE *fp_time;

  switch (build_simple_graph_instance(INPUT_FILE)) {
  case FALSE:
    printf("Input file error\n");
    return FALSE;
  case TRUE:

    // printf("init\n");

    init();

    // printf("search_by_up\n");

    result = search_by_up();
    if (result == TRUE) {

      // printf("verify\n");

      if (verify_solution() == TRUE) {
        printf(">>data: lb = %5d | ub = %5d | ", lower_bound, NB_VALUE);
        // print_solution(filename(INPUT_FILE), NB_VALUE);
      } else
        printf("Solution wrong\n");
    } else
      printf(">>data: lb = %5d | ub = %5d | ", NB_VALUE + 1, upper_bound);
    break;
  }

  // printf("ok\n");

  // a_tms = ( struct tms *) malloc( sizeof (struct tms));
  mess = times(a_tms);
  endtime = a_tms->tms_utime;
  // free(a_tms);

  printf(" NB_UNIT = %10d | conflicts = %10d | bactracks = %10d |", NB_UNIT,
         NB_BACK, NB_BRANCHE);

  printf(" time = %ld\n", (endtime - begintime));

  // fp_time = fopen("result", "a");
  // fprintf(fp_time, "color6 %s %5.3f %d %d %d %d %d %d %d %d\n",
  // 	  INPUT_FILE, ((double)(endtime-begintime)/CLK_TCK),
  // 	  NB_BACK, NB_BRANCHE, NB_UNIT, NB_NODE, NB_EDGE,
  // 	  NODE_STACK_fill_pointer==NB_NODE, NB_VALUE, H);
  // printf("color6 %s %5.3f %d %d %d %d %d %d %d %d\n",
  // 	  INPUT_FILE, ((double)(endtime-begintime)/CLK_TCK),
  // 	 NB_BACK, NB_BRANCHE, NB_UNIT, NB_NODE, NB_EDGE,
  // 	 NODE_STACK_fill_pointer==NB_NODE, NB_VALUE, H);
  // fclose(fp_time);

  // free_graph_instance();

  return result;
}

int main(int argc, char *argv[]) {
	
  long begintime, endtime;
  long mess;
  // struct tms *a_tms;

  parse_parameters(argc, argv);
  if (HELP_FLAG == TRUE) {
    printf("using the following parameters (the order does not matter)\n\n");
    printf("your input file\n");
    // printf("-nbColors N: N is the number of colors\n");
		printf("-s N: N is the random seed\n");
    printf("-l N: N is the lower bound on the number of colors\n");
    printf("-u N: N is the upper bound on the number of colors\n");
    printf("-f N : N is the format used, N=1, DIMACS format\n");
    printf("                             N=2, format simplified by Corinne\n");
    printf("-help: the current help message\n");
    return FALSE;
  }

  printf(">>statistics: lb ub NB_UNIT conflicts bactracks time\n");

  // NB_VALUE = upper_bound;

  a_tms = (struct tms *)malloc(sizeof(struct tms));
  mess = times(a_tms);
  begintime = a_tms->tms_utime;

  // if(lower_bound == upper_bound) {
  printf(">>data: lb = %5d | ub = %5d | ", lower_bound, upper_bound);
  mess = times(a_tms);
  endtime = a_tms->tms_utime;
  printf(" NB_UNIT = %10d | conflicts = %10d | bactracks = %10d |", NB_UNIT,
         NB_BACK, NB_BRANCHE);
  printf(" time = %ld\n", (endtime - begintime));
  // }

  NB_VALUE = (policy == TOP_DOWN
                  ? upper_bound - 1
                  : (policy == BOTTOM_UP ? lower_bound
                                         : (lower_bound + upper_bound) / 2));
  while (lower_bound < upper_bound) {

    // printf("c try with #colors = %d\n", NB_VALUE);

    int success = solve(begintime);

    if (success) {
      upper_bound = NB_VALUE--;

      // printf( "new UB %8d search time = %10f parse time = %10f encode time =
      // %10f conflicts = %8ld\n",
      // 				upper_bound, search_time, parse_time,
      // encode_time, n_conflicts);

    } else {
      lower_bound = ++NB_VALUE;

      // printf( "new LB %8d search time = %10f parse time = %10f encode time =
      // %10f conflicts = %8ld\n",
      // 				lower_bound, search_time, parse_time,
      // encode_time, n_conflicts);
    }

    if (policy == BINARY)
      NB_VALUE = (lower_bound + upper_bound) / 2;
  }

  free(a_tms);

  return TRUE;
}
