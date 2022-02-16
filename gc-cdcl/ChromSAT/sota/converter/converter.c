#include "converter.h"
#include <assert.h>

static Edge * new_edge(int a, int b){
	Edge * myEdge;
	int tmp;
	if(a == b){
		printf("Error: found an edge %d,%d, left must unequal to right\n", a,b);
		exit(0);
	}
	if(a < 1 || b < 1){
		printf("Error: found an edge %d,%d, with left or right below 1\n", a,b);
		exit(0);
	}
	if(a > b){
		tmp = a;
		a = b;
		b = tmp;
	}
	myEdge = malloc(sizeof(Edge));
	myEdge->a = a;
	myEdge->b = b;
	return myEdge;
}

static int calculate_variable_encoding(int node, int colour){
	return (node-1) * nr_of_colours + colour + 1;
}

void init_vars(int nrnodes, int nredges, int nrcolours, FILE * edgeFile){
	int nr_read, i,j, left, right;
	char buf[1000];
	i=0;
	nr_of_discarded_nodes = 0;
	nr_of_nodes = nrnodes;
	nr_of_edges = nredges;
	nr_of_colours = nrcolours;
	edges = malloc(nr_of_edges * sizeof(Edge));
	degree_info = malloc(nr_of_nodes * sizeof(int));
	
	for (i=0;i<nr_of_nodes;i++){
		for (j=0;j<nr_of_nodes;j++){ 
			adj[i][j] = FALSE;}}
	
	for (i=0;i<nr_of_nodes;i++) valid[i] = TRUE;
	best_clique = 0;
	num_prob = 0;
	max_prob = 10000;	
	
	for(i = 0; i < nr_of_nodes; i++){
		degree_info[i] =0;}
	i=0;	
	while (!feof(edgeFile)){
		nr_read= fscanf(daInput, "e %d %d\n", &left, &right);
		if(nr_read == 2){
			edges[i++] = *new_edge(left, right);
			neighbours[left-1][degree_info[left-1]]   = right - 1;
			neighbours[right-1][degree_info[right-1]] = left - 1;
			degree_info[left-1]++; 
			degree_info[right-1]++;
	    	adj[left-1][right-1] = TRUE;
	    	adj[right-1][left-1] = TRUE;}
		else{
			nr_read = fscanf(daInput, "c %s", buf);
			fgets( buf, sizeof(buf), daInput );
			/*
			if(nr_read != 1){

				printf("nr_read=%d: %s\n", nr_read, buf);

				printf("Parse error while readingin edges, succefully read:%d edges\n",
				        i);
				exit(0);}
			*/
		}
	}
	if(i != nr_of_edges){
		printf("Error: found %d edges expected %d edges\n", i , nr_of_edges);
		exit(1);}

	for (i =0; i < nr_of_nodes; i++) {
	  if (degree_info[i] < nr_of_colours && degree_info[i] > 0) {
	    nr_of_discarded_nodes++;
	  }
	}
	printf("We could discard %d nodes due to lack of sufficient edges..\n",
	       nr_of_discarded_nodes);
}

void free_vars(){
	free(edges);
}

void print_header(){	
	int nr_var, nr_clause, nr_atmost_one, nr_at_least_one, nr_conflict, nr_e;		
	nr_atmost_one = (nr_of_colours * (nr_of_colours-1) / 2) * nr_of_nodes;
	nr_conflict = (nr_of_edges * nr_of_colours);
	nr_at_least_one = (nr_of_nodes); 	
	nr_var = nr_of_nodes * nr_of_colours;
	nr_e = (nr_of_nodes * (nr_of_nodes-1) / 2) + 1;
	nr_var = nr_var + nr_e;
	nr_clause = nr_at_least_one + nr_atmost_one + nr_conflict + nr_of_edges + 1 + nr_of_colours;
	if(p_clique == TRUE){
		nr_clause += lb;}
	if(p_mergeliterals == TRUE){		
		fprintf(daOutput, "e info %d %d\n", nr_of_nodes, nr_of_colours);}
	fprintf(daOutput, "p cnf %d %d\n", nr_var, nr_clause);
}

void print_at_least_one(){
	int i,j;
	for(i=1;i<=nr_of_nodes;i++){
		for(j=0;j<nr_of_colours;j++){
			fprintf(daOutput, "%d ", calculate_variable_encoding(i,j));
		}
		fprintf(daOutput, "0\n");
	}
}

void print_at_most_one(){
	int i,j,k;
	for(i=1;i<=nr_of_nodes;i++){
		for(j = 1; j < nr_of_colours; j++){
			for(k=0; k < j; k++){
				fprintf(daOutput, "-%d -%d 0\n", calculate_variable_encoding(i,k), calculate_variable_encoding(i,j));
			}
		}
	}
}

void print_conflict(){
	int i,j;
	for(i=0; i < nr_of_edges;i++){
	  if (degree_info[edges[i].a -1] < nr_of_colours
	      || degree_info[edges[i].b -1] < nr_of_colours ) {
	    continue;
	  }
		for(j=0;j<nr_of_colours;j++){
			fprintf(daOutput, "-%d -%d 0\n", calculate_variable_encoding(edges[i].a, j), calculate_variable_encoding(edges[i].b, j));
		}
	}
}

void print_e(){
	int * forbidden, i , base, current_bound;
	base = calculate_variable_encoding(nr_of_nodes, nr_of_colours);
	
	fprintf(daOutput, "%d 0\n", base + (nr_of_nodes * (nr_of_nodes-1) / 2));
	
	/* ff alle aan elkaar geknoopte dingen wegstrepen */
	forbidden = calloc(nr_of_nodes * (nr_of_nodes-1) / 2, sizeof(int));
	for(i=0;i<nr_of_edges;i++){		
		forbidden[( (edges[i].b - 1) * (edges[i].b - 2) / 2 ) +  (edges[i].a - 1)] = 1;}
		
	/* nu voor iedere niet verbonden graaf een variable aan maken */
	current_bound = 2;
	for(i=0; i <  (nr_of_nodes * (nr_of_nodes-1) / 2) ; i++){
		if((current_bound) * (current_bound - 1) / 2 == i){
			current_bound++;}
		if(forbidden[i]){ 
			fprintf(daOutput, "-%d 0\n", base + i);}}
}

/* The following code is coppied from micheal trick's implemantation of the DSATUR algorithm 
   mat.gsia.cmu.edu/COLOR/solvers/trick.c */
int greedy_clique(int *valid, int *clique)
{
  int i,j,k;
  int max;
  int place,done;
  int *order;
  int weight[MAX_NODE];
  
  for (i=0;i<nr_of_nodes;i++) clique[i] = 0;
  order = (int *)calloc(nr_of_nodes+1,sizeof(int));
  place = 0;
  for (i=0;i<nr_of_nodes;i++) {
    if (valid[i]) {
      order[place] = i;
      place++;
    }
  }
  for (i=0;i<nr_of_nodes;i++)
    weight[i] = 0;
  for (i=0;i<nr_of_nodes;i++) 
    {
      if (!valid[i]) continue;
      for (j=0;j<nr_of_nodes;j++) 
	{
	  if (!valid[j]) continue;
	  if (adj[i][j]) weight[i]++;
	}
    }
  

  done = FALSE;
  while (!done) {
    done = TRUE;
    for (i=0;i<place-1;i++) {
      j = order[i];
      k = order[i+1];
      if (weight[j] < weight[k]) {
	order[i] = k;
	order[i+1] = j;
	done = FALSE;
      }
    }
  }


  clique[order[0]] = TRUE;
  for (i=1;i<place;i++) {
    j = order[i];
    for (k=0;k<i;k++) {
      if (clique[order[k]] && !adj[j][order[k]]) break;
    }
    if (k==i) {
      clique[j] = TRUE;
    }
    else clique[j] = FALSE;
    
  }
  max = 0;
  for (i=0;i<place;i++) 
    if (clique[order[i]]) max ++;
  
  free(order); 
  return max;
}

/* The following code is coppied from micheal trick's implemantation of the DSATUR algorithm
   mat.gsia.cmu.edu/COLOR/solvers/trick.c */
int max_w_clique(int * valid, int * clique, int lower, int target){
  int start,j,k;
  int incumb,new_weight;
  int *valid1,*clique1;
  int *order;
  int *value;
  int i,place,finish,done,place1;
  int total_left;
  
  num_prob++;
  if (num_prob > max_prob) return -1;
  for (j=0;j<nr_of_nodes;j++) clique[j] = 0;
  total_left = 0;
  for (i=0;i<nr_of_nodes;i++)
    if (valid[i]) total_left ++;
  if (total_left < lower) {
    return 0.0;
  }

  order = (int *)calloc(nr_of_nodes+1,sizeof(int));
  value = (int *) calloc(nr_of_nodes,sizeof(int));
  incumb = greedy_clique(valid,clique);
  if (incumb >=target) return incumb;
  if (incumb > best_clique) {
    best_clique=incumb;
  }
  
  place = 0;
  for (i=0;i<nr_of_nodes;i++) {
    if (clique[i]) {
      order[place] = i;
      total_left --;
      place++;
    }
  }
  start = place;
  for (i=0;i<nr_of_nodes;i++) {
    if (!clique[i]&&valid[i]) {
      order[place] = i;
      place++;
    }
  }
  finish = place;
  for (place=start;place<finish;place++) {
    i = order[place];
    value[i] = 0;
    for (j=0;j<nr_of_nodes;j++) {
      if (valid[j] && adj[i][j]) value[i]++;
    }
  }

  done = FALSE;
  while (!done) {
    done = TRUE;
    for (place=start;place<finish-1;place++) {
      i = order[place];
      j = order[place+1];
      if (value[i] < value[j] ) {
	order[place] = j;
	order[place+1] = i;
	done = FALSE;
      }
    }
  }
  free(value);
  for (place=start;place<finish;place++) {
    if (incumb + total_left < lower) {
      return 0;
    }
    
    j = order[place];
    total_left --;
    
    if (clique[j]) continue;
    
    valid1 = (int *)calloc(nr_of_nodes,sizeof(int));
    clique1 = (int *)calloc(nr_of_nodes,sizeof(int));
    for (place1=0;place1 < nr_of_nodes;place1++) valid1[place1] = FALSE;
    for (place1=0;place1<place;place1++) {
      k = order[place1];
      if (valid[k] && (adj[j][k])){
	valid1[k] = TRUE;
      }
      else
	valid1[k] = FALSE;
    }
    new_weight = max_w_clique(valid1,clique1,incumb-1,target-1);
    if (new_weight+1 > incumb)  {
      incumb = new_weight+1;
      for (k=0;k<nr_of_nodes;k++) clique[k] = clique1[k];
      clique[j] = TRUE;
      if (incumb > best_clique) {
	best_clique=incumb;
      }
    }

    free(valid1);
    free(clique1);
    if (incumb >=target) break;
  }
  free(order);
  return(incumb);
}

/* The following code is coppied from micheal trick's implemantation of the DSATUR algorithm
   mat.gsia.cmu.edu/COLOR/solvers/trick.c */
void find_clique(){
	lb = max_w_clique(valid,clique,0,nr_of_nodes);
}

/* The following code is coppied from micheal trick's implemantation of the DSATUR algorithm
   mat.gsia.cmu.edu/COLOR/solvers/trick.c */
void print_clique(){	  
	int daCol = 0,i,j;
	if(lb > nr_of_colours){
		printf("warning->found clique of size:%d which is larger than nr of colors...\n", lb);}
	else{
		printf("found a clique of size:%d..\n", lb);}	
	for (i=0;i<nr_of_nodes;i++){
		if (clique[i]){
			fprintf(daOutput, "%d 0\n", (i * nr_of_colours) + 1 + (daCol++));		  
			for (j=0;j<nr_of_nodes;j++){
				if ((i!=j)&&clique[j] && (!adj[i][j])){
					printf("clique is no clique\n");
					fflush(stdout);
					assert(0);}}
		}
	}
}

void print_dummies(){
	int i;
	for(i = 0; i < nr_of_colours; i++){
		fprintf(daOutput, "%d 0\n", calculate_variable_encoding(nr_of_nodes - i, (i)));
	}
}
