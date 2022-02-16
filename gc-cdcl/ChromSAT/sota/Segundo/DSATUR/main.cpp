#include "graphe.h"
#include "dsatur_algo.h"

#include <cstring>

using namespace std;

int temps_algo = 1200;

int get_problem_name(char* pname, const char* efname)
{
   int    rval = 0;
   int    len = 0;

   const char * fname = strrchr(efname,'/');
   const char * lastdot = strrchr(efname,'.');

   if(!fname) {
      /* no slashes in efname.*/
      fname = efname;
   } else {
      fname++;
   }
   if (lastdot) {
      len = lastdot - fname + 1;
   } else {
      len = strlen(fname);
   }

	int tmp = snprintf(pname,len,"%s",fname);
   if (tmp < 0) {
      rval = 1;
   }
//   printf("Extracted problem name %s\n",pname);

   return rval;
}

int main(int argc,char**argv){

	if(argc!=4 and argc!=5){
		cout<<argv[0]<<" <instancename> <VSR> <time> <seed>\n";
		cout<<"VSR : 1 for DSAT rule\n      2 for PASS rule\n";
		return 0;
	}

	C_Graphe G;
	DSATUR_ dsat_;
	ifstream fichier;
	double densite;

	//Instance name
	char *pname,*argv1_cpy;
	pname = (char *)malloc(strlen(argv[1])*sizeof(char));
	argv1_cpy = (char *)malloc(strlen(argv[1])*sizeof(char));
	strncpy(argv1_cpy,argv[1],strlen(argv[1]));
	argv1_cpy[strlen(argv[1])-1] = '\0';
	get_problem_name(pname,argv1_cpy);

	//Reading input graph
	fichier.open(argv[1], ifstream::in);
	if (!fichier){
		cerr << "Erreur à l'ouverture du fichier instance\n";
		return 1;
	}

	//Initialisation
	G.lecture(fichier);
	string s(pname);
	G.name = s;
	densite = (2*G.nb_aretes)/(double)(G.nb_sommets*(G.nb_sommets-1));
	dsat_.G = G;
	temps_algo = atoi(argv[3]);

	std::vector<int> solution(G.nb_sommets);
	
	int seed = 12345;
	if(argc == 5)
		seed = atoi(argv[4]);

	dsat_.DSATUR_algo(
	    G, temps_algo, atoi(argv[2]), 0, G.nb_sommets, seed);

	cout << "Instance : " << pname << "\n";
	cout << "UB : " << dsat_.UB << "\n";
	cout << "time : " << dsat_.time_spent << "\n";
	cout << "nodes : " << dsat_.nombre_noeuds << "\n";


	free(pname);
	free(argv1_cpy);

	return 0;
}
