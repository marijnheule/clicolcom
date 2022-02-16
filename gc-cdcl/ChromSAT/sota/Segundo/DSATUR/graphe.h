#ifndef Graphe_H
#define Graphe_H

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

class C_Graphe{
	public:
		/************************************************************/
		/********************STRUCTURE DU GRAPHE*********************/
		/************************************************************/
		unsigned int seed;
		int nb_sommets;
		int nb_aretes;
		bool init;
		string name;

		vector < vector<int> > matrice_adjacence;
		vector < list<int> > sommets_voisins;
		vector < vector<double> > valeurs_arcs;
		vector <int> existe;

		bool *adj;
		int *sommets_voisins_bis;
		int *degre;

		//Graphe comparabilite
		int FLAG;
		vector< vector<int> > CLASS;

		/************************************************************/
		/************************CONSTRUCTEUR************************/
		/************************************************************/
		//Constructeur
		C_Graphe(){
			nb_sommets = 0;
			nb_aretes = 0;
			seed = -1;
			FLAG = 0;
			init = false;
		};
		//Destructeur
		virtual ~C_Graphe();

		/************************************************************/
		/*************************FONCTIONS**************************/
		/************************************************************/
		void lecture(istream & fic);//Fonction de lecture d'un graphe sous format .col
		void afficher_graphe();//Fonction d'affichage du graphe
		void init_graphe(int nb_som, double densite);
		void init_graphe_bis(int nb_som, double densite);
		int trans_comparabilite();
		int explorer(vector< vector<int> > &classe, int u, int v, int k);
		void dijkstra(int depart, int fin, double *valeur_chemin, list<int> &chemin);
		void complementaire(C_Graphe &Comp);
		void construire_dijkstra(C_Graphe &D);
		void graphe_oriente();
		void graphe_display(int num);
		bool chemin_bool(int i, int j, int prof);
		void supprimer_voisin(int sommet, int voisin);
		void ecrire_DIMACS(char *nomfichier);
		void ecrire_DIMACS2(char *nomfichier);
		void parcours_profondeur_rec(int *visite, int sommet);
		bool connexe();

		void copier(C_Graphe &G);

		//Graphe chordal
		void maximum_cardinality_search(vector<int> &alpha, vector<int> &alpha_1);
		bool chordal(vector<int> &alpha, vector<int> &alpha_1);

		//Graphe comparabilite
		bool comparability();
		void explore(int i, int j, int k);

		//Operations sur les graphes
		void sommet_contraction(int u, int v);
		void aretes_ajout(int u, int v);
};

int creer_graphe_GA(C_Graphe &GC, C_Graphe &G_Cornaz, bool save);
int creer_graphe_GA2(C_Graphe &GC, C_Graphe &G_Cornaz, bool save, int *indices, bool **save_simplicial);
void creer_graphe_GC(C_Graphe &G, C_Graphe &GC, int *coloration_courante, int solution_courante, int *representant);
void creer_graphe_GC2(C_Graphe &G, C_Graphe &GC, int *coloration_courante, int solution_courante, int *representant,int *matrice_adjacence);

#endif
