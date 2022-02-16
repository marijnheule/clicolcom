#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "converter.h"


int main(int argc, char *argv[]){
	
	int nr_read, nrcolour, nredges, nrnodes;
	char readBuf[1000];
	char file_string[1000];
	
	if(argc != 6){
		printf("usage ./converter graphFile [nr of colours] [clique 0=no clique 1=clique] [mergeliterals 0=normal 1=mergeliterals] [cnfFile]\n");
		return 1;}
	daInput = fopen(argv[1],"r");
	if(daInput ==0){
		printf("Error: problems opening file %s...\n", argv[1]);
		return 1;}	
	nrcolour = atoi(argv[2]);
	if(nrcolour < 3){
		printf("Error: nr of colours should be greater than 2...\n");
		return 1;}
	p_clique = atoi(argv[3]);
	if (!(p_clique == TRUE || p_clique == FALSE)){
		printf("no valid print clique value provided -> clique 0=no clique 1=clique\n");
		return 1;}	
	p_mergeliterals = atoi(argv[4]);
	if (!(p_mergeliterals == TRUE || p_mergeliterals == FALSE)){
		printf("no valid print evar value provided -> evars 0=normal 1=evar\n");
		return 1;}
	sprintf(file_string, "%s.cnf", argv[5]);
	daOutput = fopen(file_string,"w");
	if(daOutput ==0){
		printf("Error: problems opening file %s.cnf...\n", argv[6]);
	return 1;}		
	
	nr_read= fscanf(daInput, "p %s %d %d\n", readBuf, &nrnodes, &nredges);

	while(nr_read != 3 && nr_read != EOF){ 
		fscanf(daInput, "%s\n", readBuf);
		
		nr_read= fscanf(daInput, "p %s %d %d\n", readBuf, &nrnodes, &nredges);
		
		readBuf[4] = '\0';
		if( strcmp(readBuf, "edge" ) )
			nr_read = 0;
	
	}	
	if(nr_read != 3){
		printf("Error: no p line found exiting...\n");
		return 1;}	
	
	
	init_vars(nrnodes+nrcolour, nredges, nrcolour, daInput);
	find_clique();		
	print_header();
	print_at_least_one();
	print_at_most_one();
	print_conflict();
	print_e();
	print_dummies();
	if(p_clique == TRUE){
		print_clique(); /*Must be the last added clauses because minisat removes satisfied clauses*/
	}
	free_vars();
	fclose(daInput);
	fclose(daOutput);
	return 0;
}
