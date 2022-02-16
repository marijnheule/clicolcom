#include "bornes.h"

//#define debug_GA
//#define write_cuts

//Variables globales pour le callback continu
double time_limit_borne = 3600;

boolean temps(int level,int i, int n,int max,double user_time ,double system_time, clique_options *){
	return (user_time < time_limit_borne);
}

int clique_cliquer_init(C_Graphe &GC, double cutoff, int *solution, int *size_solution){
	#ifdef DEBUG
	cout<<"DEBUT CLIQUER";
	#endif
	graph_t *g;
	g=(graph_t *)calloc(1,sizeof(graph_t));
	g->n = GC.nb_sommets;
	g->edges=(set_t *)calloc(g->n,sizeof(set_t));
	for (int i=0; i<g->n; i++){
		g->edges[i]=set_new(g->n);
	}
	g->weights=(int *)calloc(g->n,sizeof(int));
	for(int i=0; i<g->n; i++){
		g->weights[i]=1;
	}
	for(int i=0; i<g->n; i++){
		for(int j=i+1; j<g->n; j++){
			if(GC.adj[i*GC.nb_sommets+j] && i!=j){
				GRAPH_ADD_EDGE(g,i,j);
			}
		}
	}

//	ASSERT(graph_test(g,stderr));
	set_t s;

	time_limit_borne = 5;

	clique_options *opts;
	opts=(clique_options *)malloc(sizeof(clique_options));
	opts->time_function = temps;
	opts->output = stderr;
	opts->reorder_function = reorder_by_greedy_coloring;
	opts->reorder_map = NULL;
	opts->user_function = NULL;
	opts->user_data = NULL;
	opts->clique_list = NULL;
	opts->clique_list_length = 0;

	s = clique_find_single(g,cutoff,cutoff,FALSE,opts);
	if(s == NULL){return 0;}
//	set_print(s);
	int res = set_size(s);

	(*size_solution) = 0;
	for(int i=0 ; i<GC.nb_sommets ; i++){
		if(SET_CONTAINS(s,i)){
			solution[(*size_solution)] = i;
			(*size_solution)++;
		}
	}

	free(opts);
	set_free(s);
	graph_free(g);

	#ifdef DEBUG
	cout<<"FIN CLIQUER";
	#endif

	return res;

}
