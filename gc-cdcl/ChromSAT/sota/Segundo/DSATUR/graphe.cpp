#include "graphe.h"
#include <climits>
//#define DEBUG_dijkstra
//#define DEBUG_chordal
//#define DEBUG_auxiliaire

int index(int row, int column, int vertexNumber){return row*vertexNumber + column;}
int index_magic(int row, int column, int vertexNumber_transformedGraph){return row*vertexNumber_transformedGraph + column;}

void C_Graphe::lecture(istream & fic){
	if (!fic){
		cout<<"Erreur fichier lecture graphe"<<endl;
	}else{
		int u,v;
		string m1,m2;

		fic>>m1;
		fic>>m2;
		
		//On lit jusqu'à la description du graphe en passant les infos du début
		while (((m1!="p")&&(m2!="edge"))||((m1!="p")&&(m2!="col"))){
			m1=m2;
			fic>>m2;
		}
		//On récupère les tailles
		fic>>nb_sommets;
		fic>>nb_aretes;
		
		matrice_adjacence.resize(nb_sommets);
		sommets_voisins.resize(nb_sommets);

		adj = (bool *)malloc(nb_sommets*nb_sommets*sizeof(bool));
		sommets_voisins_bis = (int *)malloc(nb_sommets*nb_sommets*sizeof(int ));
		degre = (int *)malloc(nb_sommets*sizeof(int));

		for(int i=0 ; i<nb_sommets ; i++){
			matrice_adjacence[i].resize(nb_sommets);
			sommets_voisins[i].clear();
			degre[i] = 0;
		}
		for(int i=0 ; i<nb_sommets ; i++){
			for(int j=0 ; j<nb_sommets ; j++){
				matrice_adjacence[i][j] = 0;
				sommets_voisins_bis[i*nb_sommets+j] = 0;
				adj[i*nb_sommets+j] = false;
			}
		}

		//On crée et on rentre les arêtes dans les sommets
		int nb_a = 0;
		for (int i=0; i<nb_aretes; i++){
			fic>>m1;
			fic>>u;
			fic>>v;
			u = u-1;
			v = v-1;
			
			if(!fic) {
				nb_aretes = i;
				break;
			}
			
			if(matrice_adjacence[u][v] == 0 ){
				nb_a++;
				matrice_adjacence[u][v] = 1;
				matrice_adjacence[v][u] = 1;

				adj[u*nb_sommets+v] = true;
				adj[v*nb_sommets+u] = true;

				sommets_voisins[u].push_back(v);
				sommets_voisins[v].push_back(u);

				sommets_voisins_bis[(u*nb_sommets)+degre[u]] = v;
				degre[u]++;
				sommets_voisins_bis[(v*nb_sommets)+degre[v]] = u;
				degre[v]++;
			}
//			else{cout<<"doublon\n";}
		}

//		for(int i=0 ; i<nb_sommets ; i++){
//			int j=0;
//			for(list<int>::iterator it=sommets_voisins[i].begin() ; it!=sommets_voisins[i].end() ; it++){
//				cout<<*it<<"|";
//				cout<<sommets_voisins_bis[(i*nb_sommets)+j]<<" ";
//				if(*it != sommets_voisins_bis[(i*nb_sommets)+j]){cout<<"probleme\n";}
//				j++;
//			}
//			cout<<"\n";
//			cin.get();
//		}
		for(int i=0 ; i<nb_sommets ; i++){
			for(int j=0 ; j<nb_sommets ; j++){
if(matrice_adjacence[i][j] == 1 && !adj[i*nb_sommets+j]){cout<<"problème\n";}
if(matrice_adjacence[i][j] == 0 && adj[i*nb_sommets+j]){cout<<"problème\n";}
			}
		}
		nb_aretes = nb_a;
		existe.resize(nb_sommets);
		for(int i=0 ; i<nb_sommets ; i++){
			existe[i] = 1;
		}
		for(int i=0 ; i<nb_sommets ; i++){
			sommets_voisins[i].sort();
		}
	}
	init = true;

	return;
}

void C_Graphe::parcours_profondeur_rec(int *visite, int sommet){
	if(visite[sommet] == 1){
		return;
	}else{
		visite[sommet] = 1;
		for(list<int>::iterator it_voisin = sommets_voisins[sommet].begin() ; it_voisin != sommets_voisins[sommet].end() ; it_voisin++){
			parcours_profondeur_rec(visite,*it_voisin);
		}
	}
	return;
}

bool C_Graphe::connexe(){
	int *visite = (int*)malloc(nb_sommets*sizeof(int));
	for(int i=0 ; i<nb_sommets ; i++){
		visite[i] = 0;
	}

	parcours_profondeur_rec(visite,0);

	for(int i=0 ; i<nb_sommets ; i++){
		if(visite[i] == 0){
			free(visite);
			return false;
		}
	}

	free(visite);
	return true;
}

void C_Graphe::init_graphe_bis(int nb_som, double densite){
	unsigned int s = time(NULL);
	srand(s);
	seed = s;
	nb_sommets = nb_som;
	matrice_adjacence.resize(nb_sommets);
	sommets_voisins.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		matrice_adjacence[i].resize(nb_sommets);
		sommets_voisins[i].clear();
	}
	int m = densite*nb_sommets*(nb_sommets-1)*0.5;
	int random_1,random_2;
	int cpt = 0;
	for(int i=0;i<m;i++){
		bool _OK= false;
		do{
			random_1=(int)rand();
			random_1=(random_1%(nb_sommets-1));

			random_2=(int)rand();
			random_2=(random_2%(nb_sommets-1));

			if(random_1!=random_2 && matrice_adjacence[random_1][random_2] == 0){
				matrice_adjacence[random_1][random_2] = 1;
				matrice_adjacence[random_2][random_1] = 1;
				sommets_voisins[random_1].push_back(random_2);
				sommets_voisins[random_2].push_back(random_1);
				nb_aretes++;
				_OK=true;
			}
			cpt++;
		}while(_OK==false && cpt <= 100);
	}
	return;
}

void C_Graphe::init_graphe(int nb_som, double densite){
	srand(time(NULL));
	nb_sommets = nb_som;
	matrice_adjacence.resize(nb_sommets);
	sommets_voisins.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		matrice_adjacence[i].resize(nb_sommets);
		sommets_voisins[i].clear();
	}
	double rand_num;
	for(int i=0 ; i<nb_sommets-1 ; i++){
		matrice_adjacence[i][i] = 0;
		for(int j=i+1 ; j<nb_sommets ; j++){
			rand_num = (double)rand() / (double)RAND_MAX;
			if(rand_num <= densite){
				matrice_adjacence[i][j] = 1;
				matrice_adjacence[j][i] = 1;
				sommets_voisins[j].push_back(i);
				sommets_voisins[i].push_back(j);
				nb_aretes++;
			}else{
				matrice_adjacence[i][j] = 0;
				matrice_adjacence[j][i] = 0;
			}
		}
	}
	matrice_adjacence[nb_sommets-1][nb_sommets-1] = 0;
	return;
}

void C_Graphe::afficher_graphe(){
	cout<<"Le graphe a "<<nb_sommets<<" sommets et "<<nb_aretes<<" aretes\n";
	int cpt_arete = 0;
	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 1){
				cout<<"Arete : ("<<i<<","<<j<<")\n";
				cpt_arete++;
			}
		}
	}
	cout<<(cpt_arete == nb_aretes?"Compte bon":"Erreur aretes")<<"\n";
	return;
}

int C_Graphe::trans_comparabilite(){
	int cpt = 0;	
	vector< vector<int> > classe;
	classe.resize(nb_sommets);
	int undefined = nb_sommets*2;
	for(int i=0 ; i<nb_sommets ; i++){
		classe[i].resize(nb_sommets);
		for(int j=0 ; j<nb_sommets ; j++){
			classe[i][j] = 0;
			if(matrice_adjacence[i][j] == 1){
				classe[i][j] = undefined;
			}
		}
	}
	int k=0;
	for(int u=0 ; u<nb_sommets-1 ; u++){
		for(int v=u+1 ; v<nb_sommets ; v++){if(matrice_adjacence[u][v] == 1){
			if(classe[u][v] == undefined){
				k++;
				classe[u][v] = k;
				classe[v][u] = -k;
				cpt += explorer(classe,u,v,k);
			}
		}}
	}
	return cpt;
}

int C_Graphe::explorer(vector< vector<int> > &classe, int u, int v, int k){
	int cpt = 0;
	int undefined = nb_sommets*2;
	for(list<int>::iterator it_vois = sommets_voisins[u].begin() ; it_vois != sommets_voisins[u].end() ; it_vois++){if(matrice_adjacence[*it_vois][v] == 0 || abs(classe[v][*it_vois]) < k){
		if(classe[u][*it_vois] == undefined){
			classe[u][*it_vois] = k;
			classe[*it_vois][u] = -k;
			cpt += explorer(classe,u,*it_vois,k);
		}else{
			if(classe[u][*it_vois] == -k){
				classe[u][*it_vois] = k;
				cout<<"Enlève ("<<v<<","<<*it_vois<<")\n";
				supprimer_voisin(v,*it_vois);
				supprimer_voisin(*it_vois,v);
				it_vois = sommets_voisins[u].begin();
				cpt +=explorer(classe,u,*it_vois,k);
			}
		}
	}}
	for(list<int>::iterator it_vois = sommets_voisins[v].begin() ; it_vois != sommets_voisins[v].end() ; it_vois++){if(matrice_adjacence[*it_vois][u] == 0 || abs(classe[u][*it_vois]) < k){
		if(classe[*it_vois][v] == undefined){
			classe[*it_vois][v] = k;
			classe[v][*it_vois] = -k;
			cpt += explorer(classe,*it_vois,u,k);
		}else{
			if(classe[*it_vois][v] == -k){
				classe[*it_vois][v] = k;
				cout<<"Enlève ("<<u<<","<<*it_vois<<")\n";
				supprimer_voisin(u,*it_vois);
				supprimer_voisin(*it_vois,u);
				it_vois = sommets_voisins[v].begin(); 
				cpt +=explorer(classe,*it_vois,v,k);
			}
		}
	}}
	return cpt;
}

void C_Graphe::dijkstra(int depart, int fin, double *valeur_chemin, list<int> &chemin){
	vector<int> marque;
	vector<double> valeur;
	vector<int> predecesseur;
	int compteur_sommets,choix_som;
	double distance,min;

//	cout<<"Dijkstra pcc entre "<<depart<<" et "<<fin<<"\n";

	if(sommets_voisins[depart].size() == 0 || sommets_voisins[fin].size() == 0){
		*valeur_chemin = nb_sommets;
		chemin.clear();
		return;
	}

	marque.resize(nb_sommets);
	valeur.resize(nb_sommets);
	predecesseur.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		marque[i] = 0;
		valeur[i] = nb_sommets;
		predecesseur[i] = -1;
	}
	compteur_sommets = 0;
	valeur[depart] = 0;
	while(compteur_sommets < nb_sommets){
		//choix du sommet
		choix_som=0;
		min = nb_sommets+1;
		for(int i=0 ; i<nb_sommets ; i++){
			if(marque[i] == 0 && valeur[i] < min){
				choix_som = i;
				min = valeur[i];
			}
		}
		#ifdef DEBUG_dijkstra
		cout<<"Choix : "<<choix_som<<" "<<valeur[choix_som]<<"\n";
		#endif
		//on marque le nouveau sommet considéré
		marque[choix_som] = 1;
		compteur_sommets++;
		for(list<int>::iterator it_voisin = sommets_voisins[choix_som].begin() ; it_voisin != sommets_voisins[choix_som].end() ; it_voisin++){if(marque[*it_voisin] == 0 && *it_voisin != predecesseur[*it_voisin] != choix_som){
			distance = valeur[choix_som]+valeurs_arcs[*it_voisin][choix_som];
			#ifdef DEBUG_dijkstra
			cout<<"Valeur candidate : "<<distance<<"\n";
			#endif
			if(distance < valeur[*it_voisin]){
				predecesseur[*it_voisin] = choix_som;
				valeur[*it_voisin] = distance;
				#ifdef DEBUG_dijkstra
				cout<<"Update : "<<distance<<" sommet "<<*it_voisin<<" predecesseur : "<<choix_som<<"\n";
				#endif
			}
		}
	}}
	#ifdef DEBUG_dijkstra
	cout<<"Données après algo :\n";
	for(int i=0 ; i<nb_sommets ; i++){
		cout<<"Sommet "<<i<<" "<<marque[i]<<" "<<valeur[i]<<" "<<predecesseur[i]<<"\n";
	}
	#endif
//	if (fin == 44 || depart == 44){
//	cout<<"Données après algo :\n";
//	for(int i=0 ; i<nb_sommets ; i++){
//		cout<<"Sommet "<<i<<" "<<marque[i]<<" "<<valeur[i]<<" "<<predecesseur[i]<<" "<<valeurs_arcs[i][predecesseur[i]]<<"\n";
//	}
//	}


//	reconstruction du chemin
	chemin.clear();
	int tmp_predecesseur = predecesseur[fin];
	chemin.push_front(fin);
	while(tmp_predecesseur != depart){
		if(tmp_predecesseur == -1){
			*valeur_chemin = nb_sommets;
			chemin.clear();
			return;
		}
		chemin.push_front(tmp_predecesseur);
		tmp_predecesseur = predecesseur[tmp_predecesseur];
		#ifdef DEBUG_dijkstra
		cout<<tmp_predecesseur<<"\n";
		#endif
	}
	*valeur_chemin = valeur[fin];
	chemin.push_front(depart);
//	cout<<"Fin Dijkstra\n";
	return;
}

void C_Graphe::construire_dijkstra(C_Graphe &G_dijkstra){
	G_dijkstra.nb_sommets = nb_sommets*2;
	G_dijkstra.nb_aretes = 0;
	G_dijkstra.sommets_voisins.resize(nb_sommets*2);
	G_dijkstra.valeurs_arcs.resize(nb_sommets*2);
	G_dijkstra.matrice_adjacence.resize(nb_sommets*2);

	for(int i=0 ; i<(int)G_dijkstra.valeurs_arcs.size() ; i++){
		G_dijkstra.valeurs_arcs[i].resize(G_dijkstra.nb_sommets*2);
		G_dijkstra.matrice_adjacence[i].resize(G_dijkstra.nb_sommets*2);
	}
	for(int i=0; i<G_dijkstra.matrice_adjacence.size() - 1; i++) {
		for(int j=i+1 ; j<G_dijkstra.matrice_adjacence.size() ; j++){
			G_dijkstra.matrice_adjacence[i][j] = 0;
			G_dijkstra.matrice_adjacence[j][i] = 0;
			G_dijkstra.valeurs_arcs[i][j] = 0;
			G_dijkstra.valeurs_arcs[j][i] = 0;
		}
	}
	for(int i=0; i<matrice_adjacence.size(); i++) {
		for(int j=0 ; j<matrice_adjacence.size() ; j++){
			if(matrice_adjacence[i][j] == 1){
				G_dijkstra.matrice_adjacence[i][nb_sommets+j] = 1;
				G_dijkstra.matrice_adjacence[nb_sommets+i][j] = 1;

				G_dijkstra.sommets_voisins[i].push_back(nb_sommets+j);
				G_dijkstra.sommets_voisins[nb_sommets+i].push_back(j);
//				cout<<"ajout arete ("<<i<<","<<nb_sommets+j<<") "<<"("<<nb_sommets+i<<","<<j<<")\n";

//				G_dijkstra.valeurs_arcs[i][nb_sommets+j] = 1 - y[i] - y[j];
//				G_dijkstra.valeurs_arcs[nb_sommets+i][j] = 1 - y[i] - y[j];

				G_dijkstra.nb_aretes++;
			}
		}
	}
	return ;
}

void C_Graphe::complementaire(C_Graphe &C){
	C.nb_sommets = nb_sommets;
	C.nb_aretes = 0;

	C.sommets_voisins.resize(C.nb_sommets);
	C.matrice_adjacence.resize(C.nb_sommets);
	for(int i=0 ; i<C.nb_sommets ; i++){
		C.matrice_adjacence[i].resize(nb_sommets);
	}
	
	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 0){
				C.matrice_adjacence[i][j] = 1;
				C.matrice_adjacence[j][i] = 1;
				C.nb_aretes++;
				C.sommets_voisins[i].push_back(j);
				C.sommets_voisins[j].push_back(i);
			}
		}
	}

	return;
}

bool C_Graphe::chemin_bool(int i, int j, int prof){
	if(prof > nb_sommets){
		return false;
	}
	if(i == j){
		return true;
	}
	for(list<int>::iterator it_voisin = sommets_voisins[i].begin() ; it_voisin != sommets_voisins[i].end() ; it_voisin++){
		if(chemin_bool(*it_voisin,j,prof+1)){
			return true;
		}
	}
	return false;
}

void C_Graphe::graphe_oriente(){
	ofstream fichier_graphe;
	fichier_graphe.open("digraph.dot", ios::out);

	fichier_graphe<<"digraph G {\n";
	for(int i=0 ; i<nb_sommets ; i++){
		for(int j=0 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 1){
				fichier_graphe<<"\""<<i<<"\" -> \""<<j<<"\"\n";
			}
		}
	}
	fichier_graphe<<"}\n";
	fichier_graphe.close();
	if(system("dot -Tps -o digraph.ps digraph.dot")){
		cout<<"Problème dessin\n";
	}
	return;
}

void C_Graphe::graphe_display(int num){
	ofstream fichier_graphe;

	ostringstream namebuff_dot,namebuff_ps,commande;
	namebuff_dot.str("");namebuff_dot<<"graphe_"<<num<<".dot";
	namebuff_ps.str("");namebuff_ps<<"graphe_"<<num<<".ps";

	fichier_graphe.open(namebuff_dot.str().c_str(), ios::out);

	fichier_graphe<<"graph G {\n";
	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 1){
				fichier_graphe<<"\""<<i<<"\" -- \""<<j<<"\"\n";
			}
		}
	}
	fichier_graphe<<"}\n";
	fichier_graphe.close();

	commande.str("");commande<<"dot -Tps -o "<<namebuff_ps.str().c_str()<<" "<<namebuff_dot.str().c_str();
	if(system(commande.str().c_str())){
		cout<<"Problème dessin\n";
	}
	return;
}

void C_Graphe::supprimer_voisin(int sommet, int voisin){
	matrice_adjacence[sommet][voisin] = 0;
	for(list<int>::iterator it_vois = sommets_voisins[sommet].begin() ; it_vois != sommets_voisins[sommet].end() ; it_vois++){
		if(*it_vois == voisin){
			sommets_voisins[sommet].erase(it_vois);
			return ;
		}
	}
	return ;
}

void C_Graphe::ecrire_DIMACS(char *nomfichier){
	ofstream fichier_graphe;
	fichier_graphe.open(nomfichier, ios::out);
	fichier_graphe<<"c "<<nomfichier<<"\n";
	fichier_graphe<<"c Random Graph d="<<(2*nb_aretes)/(double)(nb_sommets*(nb_sommets-1))<<" seed="<<seed<<"\n";
	fichier_graphe<<"p edge "<<nb_sommets<<" "<<nb_aretes<<"\n";
	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 1){
				fichier_graphe<<"e "<<min(i,j)+1<<" "<<max(i,j)+1<<"\n";
			}
		}
	}
	fichier_graphe.close();
	return;
}

void C_Graphe::ecrire_DIMACS2(char *nomfichier){
	ofstream fichier_graphe;
	fichier_graphe.open(nomfichier, ios::out);
	fichier_graphe<<nb_sommets<<" "<<nb_aretes<<"\n";
	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if(matrice_adjacence[i][j] == 1){
				fichier_graphe<<min(i,j)+1<<" "<<max(i,j)+1<<"\n";
			}
		}
	}
	fichier_graphe.close();
	return;
}

int ind(int row, int column, int vertexNumber){return row*vertexNumber + column;}
int ind_magic(int row, int column, int vertexNumber_transformedGraph){return row*vertexNumber_transformedGraph + column;}

int creer_graphe_GA(C_Graphe &GC, C_Graphe &G_Cornaz, bool save){
	#ifdef DEBUG_auxiliaire
	cout<<"DEBUT creer_graphe_GA\n";
	#endif
	G_Cornaz.init = true;

	int i,j,k;

	int vertexNumber = GC.nb_sommets;
	long long int vertexNumber_transformedGraph = ((GC.nb_sommets*(GC.nb_sommets-1))/2) - GC.nb_aretes;
	if(vertexNumber_transformedGraph*vertexNumber_transformedGraph > INT_MAX){
		G_Cornaz.init = false;
	#ifdef DEBUG_auxiliaire
	cout<<"Graphe trop grand\n";
	cout<<"FIN creer_graphe_GA\n";
	#endif
		return 0;
	}

	#ifdef DEBUG_auxiliaire
	cout<<"nb sommets = "<<vertexNumber_transformedGraph<<" "<<INT_MAX-(vertexNumber_transformedGraph*vertexNumber_transformedGraph)<<" "<<INT_MAX<<" "<<vertexNumber_transformedGraph*vertexNumber_transformedGraph<<"\n";
	#endif

	int *from_i_j_to_node;
	from_i_j_to_node=(int*)malloc(vertexNumber*vertexNumber*sizeof(int));
	int counter=0;
	for(i=0; i<vertexNumber; i++){from_i_j_to_node[index(i,i,vertexNumber)]=-1;}
	for(i=0; i<vertexNumber; i++){
		for(j=i+1; j<vertexNumber; j++){
			if(!GC.adj[i*GC.nb_sommets+j]){
				from_i_j_to_node[index(i,j,vertexNumber)]=counter;
				from_i_j_to_node[index(j,i,vertexNumber)]=counter;
				counter++;
			}
			else{
				from_i_j_to_node[index(i,j,vertexNumber)]=-1;
				from_i_j_to_node[index(j,i,vertexNumber)]=-1;
			}
		}
	}
//	for(i=0; i<vertexNumber-1; i++){
//		for(j=i+1; j<vertexNumber; j++){
//			if(from_i_j_to_node[index(i,j,vertexNumber)] >= -1){
//				cout<<"("<<i<<","<<j<<") = "<<from_i_j_to_node[index(i,j,vertexNumber)]<<"\n";
//			}
//		}
//	}
//	if(save){
//		for(i=0; i<vertexNumber; i++){
//			for(j=0; j<vertexNumber; j++){
//				indices[index(i,j,vertexNumber)] = from_i_j_to_node[index(i,j,vertexNumber)];
//			}
//		}
////		for(i=0; i<vertexNumber-1; i++){
////			for(j=i+1; j<vertexNumber; j++){
////				if(from_i_j_to_node[index(i,j,G.nb_sommets)] != -1){
////					cout<<"Arete ("<<i<<","<<j<<") -> "<<from_i_j_to_node[index(i,j,G.nb_sommets)]<<"\n";
////				}
////			}
////		}

////		for(i=0; i<vertexNumber-1; i++){
////			for(j=i+1; j<vertexNumber; j++){
////				if(G.matrice_adjacence[i][j] == 0 && indices[index(i,j,G.nb_sommets)] == -1){
////					cout<<"Probleme\n";
////					cin.get();
////				}
////			}
////		}
//	}

	int *order=(int*)malloc(vertexNumber*sizeof(int));
	int *anti_order=(int*)malloc(vertexNumber*sizeof(int));
	for(i=0; i<vertexNumber; i++){order[i]=i;anti_order[order[i]]=i;}

	#ifdef DEBUG_auxiliaire
	cout<<"DEBUT adjacence\n";
	#endif

	int* adjacencyMatrix;
	adjacencyMatrix=(int*)malloc(vertexNumber*vertexNumber*sizeof(int));
	for(i=0; i<vertexNumber; i++){
		for(j=0; j<vertexNumber; j++){
			adjacencyMatrix[index(i,j,vertexNumber)] = (GC.adj[i*GC.nb_sommets+j]?1:0);
		}
	}

	int* adjacencyMatrix_transformedGraph;
	adjacencyMatrix_transformedGraph=(int*)malloc(vertexNumber_transformedGraph*vertexNumber_transformedGraph*sizeof(int));
	if(adjacencyMatrix_transformedGraph == NULL){
		G_Cornaz.init = false;
		free(order);
		free(anti_order);
		free(from_i_j_to_node);
		free(adjacencyMatrix);
		free(adjacencyMatrix_transformedGraph);
	#ifdef DEBUG_auxiliaire
	cout<<"Graphe trop grand\n";
	cout<<"FIN creer_graphe_GA\n";
	#endif
		return 0;
	}

//	if(save){
//	save_simplicial.resize(vertexNumber);
//	for(i=0 ; i<vertexNumber ; i++){
//		save_simplicial[i].resize(vertexNumber);
//		for(j=0; j<vertexNumber; j++){
//			save_simplicial[i][j].resize(vertexNumber);
//			for(int k=0 ; k<vertexNumber ; k++){
//				save_simplicial[i][j][k] = 0;
//			}
//		}
//	}
//	}
	for(i=0; i<vertexNumber_transformedGraph; i++){
		for(j=0; j<vertexNumber_transformedGraph; j++){
			adjacencyMatrix_transformedGraph[index_magic(i,j,vertexNumber_transformedGraph)]=-1;
		}
	}
	int i_order;
	int j_order;
	int k_order;

	int edgeNumber_transformedGraph=0;

	for(i=0; i<vertexNumber; i++){
		i_order=order[i];
		for(j=0; j<vertexNumber; j++){
			j_order=order[j];
			if(i_order==j_order){continue;}
			for(k=0; k<vertexNumber; k++){
				k_order=order[k];
				if(i_order==k_order){continue;}
				if(j_order>=k_order){continue;}
				if(adjacencyMatrix[index(i_order,j_order,vertexNumber)]>0){continue;}//Si (i,j) n'est pas une arete
				if(adjacencyMatrix[index(i_order,k_order,vertexNumber)]>0){continue;}//Si (i,k) n'est pas une arete
				if((anti_order[i_order]<anti_order[j_order]) && (anti_order[i_order]<anti_order[k_order]) &&  (adjacencyMatrix[index(j_order,k_order,vertexNumber)]==0)){//Si la paire est simpliciale
//					if(save){
//						save_simplicial[i_order][j_order][k_order]=1;
////						save_simplicial[i_order][k_order][j_order]=1;
////						save_simplicial[j_order][i_order][k_order]=1;
////						save_simplicial[j_order][k_order][i_order]=1;
////						save_simplicial[k_order][i_order][j_order]=1;
////						save_simplicial[k_order][j_order][j_order]=1;
//					}
					continue;
				}

	//printf("%d,%d,%d--%d,%d,%d--%d,%d,%d\n",i,j,k,i_order,j_order,k_order,anti_order[i_order],anti_order[j_order],anti_order[k_order]);
	//add a new edge i -- j

				adjacencyMatrix_transformedGraph[index_magic(from_i_j_to_node[index(i_order,j_order,vertexNumber)],from_i_j_to_node[index(i_order,k_order,vertexNumber)],vertexNumber_transformedGraph)]=edgeNumber_transformedGraph;
				adjacencyMatrix_transformedGraph[index_magic(from_i_j_to_node[index(i_order,k_order,vertexNumber)],from_i_j_to_node[index(i_order,j_order,vertexNumber)],vertexNumber_transformedGraph)]=edgeNumber_transformedGraph;
				edgeNumber_transformedGraph++;
			}
		}
	}

	#ifdef DEBUG_auxiliaire
	cout<<"FIN adjacence\n";
	#endif

	G_Cornaz.nb_sommets = vertexNumber_transformedGraph;
	G_Cornaz.nb_aretes = edgeNumber_transformedGraph;
	
	G_Cornaz.adj = (bool *)malloc(vertexNumber_transformedGraph*vertexNumber_transformedGraph*sizeof(int));
	bool res;
	for(i=0; i<vertexNumber_transformedGraph-1; i++){
		G_Cornaz.adj[i*G_Cornaz.nb_sommets+i] = false;
		for(j=i+1; j<vertexNumber_transformedGraph; j++){
			res = (adjacencyMatrix_transformedGraph[index_magic(i,j,vertexNumber_transformedGraph)] >= 0);
			G_Cornaz.adj[i*G_Cornaz.nb_sommets+j] = res;
			G_Cornaz.adj[j*G_Cornaz.nb_sommets+i] = res;
		}
	}
	free(order);
	free(anti_order);
	free(from_i_j_to_node);
	free(adjacencyMatrix);
	free(adjacencyMatrix_transformedGraph);
	return 1;
}

int creer_graphe_GA2(C_Graphe &GC, C_Graphe &G_Cornaz, bool save, int *indices, bool **save_simplicial){
	#ifdef DEBUG_auxiliaire
	cout<<"DEBUT creer_graphe_GA\n";
	#endif
	G_Cornaz.init = true;

	int i,j,k;

	int vertexNumber = GC.nb_sommets;
	long long int vertexNumber_transformedGraph = ((GC.nb_sommets*(GC.nb_sommets-1))/2) - GC.nb_aretes;
	if(vertexNumber_transformedGraph*vertexNumber_transformedGraph > INT_MAX){
		G_Cornaz.init = false;
	#ifdef DEBUG_auxiliaire
	cout<<"Graphe trop grand\n";
	cout<<"FIN creer_graphe_GA\n";
	#endif
		return 0;
	}

	#ifdef DEBUG_auxiliaire
	cout<<"nb sommets = "<<vertexNumber_transformedGraph<<" "<<INT_MAX-(vertexNumber_transformedGraph*vertexNumber_transformedGraph)<<" "<<INT_MAX<<" "<<vertexNumber_transformedGraph*vertexNumber_transformedGraph<<"\n";
	#endif

	int *from_i_j_to_node;
	from_i_j_to_node=(int*)malloc(vertexNumber*vertexNumber*sizeof(int));
	int counter=0;
	for(i=0; i<vertexNumber; i++){from_i_j_to_node[index(i,i,vertexNumber)]=-1;}
	for(i=0; i<vertexNumber; i++){
		for(j=i+1; j<vertexNumber; j++){
			if(!GC.adj[i*GC.nb_sommets+j]){
				from_i_j_to_node[index(i,j,vertexNumber)]=counter;
				from_i_j_to_node[index(j,i,vertexNumber)]=counter;
				counter++;
			}
			else{
				from_i_j_to_node[index(i,j,vertexNumber)]=-1;
				from_i_j_to_node[index(j,i,vertexNumber)]=-1;
			}
		}
	}
	if(save){
		for(i=0; i<vertexNumber; i++){
			for(j=0; j<vertexNumber; j++){
				indices[index(i,j,vertexNumber)] = from_i_j_to_node[index(i,j,vertexNumber)];
			}
		}
//		for(i=0; i<vertexNumber-1; i++){
//			for(j=i+1; j<vertexNumber; j++){
//				if(!GC.adj[i*GC.nb_sommets+j] && indices[index(i,j,GC.nb_sommets)] == -1){
//					cout<<"Probleme\n";
//					cin.get();
//				}
//			}
//		}
	}

	int *order=(int*)malloc(vertexNumber*sizeof(int));
	int *anti_order=(int*)malloc(vertexNumber*sizeof(int));
	for(i=0; i<vertexNumber; i++){order[i]=i;anti_order[order[i]]=i;}

	#ifdef DEBUG_auxiliaire
	cout<<"DEBUT adjacence\n";
	#endif

	bool* adjacencyMatrix_transformedGraph;
	adjacencyMatrix_transformedGraph=(bool*)malloc(vertexNumber_transformedGraph*vertexNumber_transformedGraph*sizeof(bool));
	if(adjacencyMatrix_transformedGraph == NULL){
		G_Cornaz.init = false;
		free(order);
		free(anti_order);
		free(from_i_j_to_node);
		free(adjacencyMatrix_transformedGraph);
	#ifdef DEBUG_auxiliaire
	cout<<"Graphe trop grand\n";
	cout<<"FIN creer_graphe_GA\n";
	#endif
		return 0;
	}

	if(save){
	for(i=0 ; i<vertexNumber ; i++){
		for(j=0; j<vertexNumber; j++){
			for(k=0 ; k<vertexNumber ; k++){
				save_simplicial[(i*vertexNumber)+j][k] = false;
			}
		}
	}
	}

	for(i=0; i<vertexNumber_transformedGraph; i++){
		for(j=0; j<vertexNumber_transformedGraph; j++){
			adjacencyMatrix_transformedGraph[index_magic(i,j,vertexNumber_transformedGraph)]=false;
		}
	}
	int i_order;
	int j_order;
	int k_order;

	int edgeNumber_transformedGraph=0;

	for(i=0; i<vertexNumber; i++){
		i_order=order[i];
		for(j=0; j<vertexNumber; j++){
			j_order=order[j];
			if(i_order==j_order){continue;}
			for(k=0; k<vertexNumber; k++){
				k_order=order[k];
				if(i_order==k_order){continue;}
				if(j_order>=k_order){continue;}
				if(GC.adj[index(i_order,j_order,vertexNumber)]){continue;}//Si (i,j) n'est pas une arete
				if(GC.adj[index(i_order,k_order,vertexNumber)]){continue;}//Si (i,k) n'est pas une arete
				if((anti_order[i_order]<anti_order[j_order]) && (anti_order[i_order]<anti_order[k_order]) &&  !(GC.adj[index(j_order,k_order,vertexNumber)])){//Si la paire est simpliciale
					if(save){
						save_simplicial[(i_order*vertexNumber)+j_order][k_order]=true;
//						save_simplicial[i_order][k_order][j_order]=1;
//						save_simplicial[j_order][i_order][k_order]=1;
//						save_simplicial[j_order][k_order][i_order]=1;
//						save_simplicial[k_order][i_order][j_order]=1;
//						save_simplicial[k_order][j_order][j_order]=1;
					}
					continue;
				}
				adjacencyMatrix_transformedGraph[index_magic(from_i_j_to_node[index(i_order,j_order,vertexNumber)],from_i_j_to_node[index(i_order,k_order,vertexNumber)],vertexNumber_transformedGraph)]=true;
				adjacencyMatrix_transformedGraph[index_magic(from_i_j_to_node[index(i_order,k_order,vertexNumber)],from_i_j_to_node[index(i_order,j_order,vertexNumber)],vertexNumber_transformedGraph)]=true;
				edgeNumber_transformedGraph++;
			}
		}
	}

	#ifdef DEBUG_auxiliaire
	cout<<"FIN adjacence\n";
	#endif

	G_Cornaz.nb_sommets = vertexNumber_transformedGraph;
	G_Cornaz.nb_aretes = edgeNumber_transformedGraph;
	
	G_Cornaz.adj = (bool *)malloc(vertexNumber_transformedGraph*vertexNumber_transformedGraph*sizeof(bool));
	bool res;
	for(i=0; i<vertexNumber_transformedGraph-1; i++){
		G_Cornaz.adj[i*G_Cornaz.nb_sommets+i] = false;
		for(j=i+1; j<vertexNumber_transformedGraph; j++){
			res = (adjacencyMatrix_transformedGraph[index_magic(i,j,vertexNumber_transformedGraph)]);
			G_Cornaz.adj[i*G_Cornaz.nb_sommets+j] = res;
			G_Cornaz.adj[j*G_Cornaz.nb_sommets+i] = res;
		}
	}
	free(order);
	free(anti_order);
	free(from_i_j_to_node);
	free(adjacencyMatrix_transformedGraph);
	return 1;
}


void creer_graphe_GC(C_Graphe &G, C_Graphe &GC, int *coloration_courante, int solution_courante, int *representant){
	int cpt_non_colorie;
	cpt_non_colorie = 0;
	//On construit un sommet par groupe de couleurs pour remplacer les sommets déja coloriés en élisant un représentant
	representant[G.nb_sommets+2-1] = 1;

	//Vecteur tampon matrice_adjacence
	vector< vector<int> > matrice_adjacence;
	matrice_adjacence.resize(G.nb_sommets);
	for(int i=0 ; i<G.nb_sommets ; i++){
		matrice_adjacence[i].resize(G.nb_sommets);
		if(coloration_courante[i] == -1)
			cpt_non_colorie++;
	}

	int *representant_couleur;
	representant_couleur = (int *)malloc(sizeof(int)*(solution_courante+1));
	for(int i=0 ; i<G.nb_sommets ; i++){
		representant[i] = -1;
	}
	for(int i=1 ; i<=solution_courante ; i++){
		representant_couleur[i] = -1;
	}
	for(int i=0 ; i<G.nb_sommets ; i++){
		#ifdef DEBUG
		cout<<i<<" "<<coloration_courante[i]<<" "<<representant_couleur[coloration_courante[i]]<<"\n";
		#endif
		if(coloration_courante[i] != -1 && representant_couleur[coloration_courante[i]] == -1){
			representant_couleur[coloration_courante[i]] = i;
		}
	}
	for(int i=0 ; i<G.nb_sommets ; i++){
		if(coloration_courante[i] == -1){
			representant[i] = i;
		}else{
			representant[i] = representant_couleur[coloration_courante[i]];
		}
	}
	#ifdef DEBUG
	cout<<"Representants couleur\n";
	for(int i=1 ; i<=solution_courante ; i++){
		cout<<"Couleur "<<i<<" représentée par "<<representant_couleur[i]<<"\n";
	}
	cout<<"Representants\n";
	for(int i=0 ; i<G.nb_sommets ; i++){
		cout<<"Sommet "<<i<<" représenté par "<<representant[i]<<"\n";
	}
	#endif
	//Initialisation
	for(int i=0 ; i<G.nb_sommets ; i++){
		for(int j=i ; j<G.nb_sommets ; j++){
			matrice_adjacence[i][j] = 0;
			matrice_adjacence[j][i] = 0;
		}
	}
	//On change les vecteurs d'adjacence des représentants
	for(int i=0 ; i<G.nb_sommets ; i++){
		for(int j=i ; j<G.nb_sommets ; j++){
			if(G.adj[i*G.nb_sommets+j] == 1){
				if(matrice_adjacence[representant[i]][representant[j]] == 0){GC.nb_aretes++;}
				matrice_adjacence[representant[i]][representant[j]] = 1;
				matrice_adjacence[representant[j]][representant[i]] = 1;
			}
		}
	}
	//On crée finalement la matrice d'adjacence de GC
	vector<int> new_indice;
	new_indice.resize(G.nb_sommets);
	int cpt_indice = 0;
	for(int i=0 ; i<G.nb_sommets ; i++){
		if(representant[i] == i){
			new_indice[i] = cpt_indice;
			cpt_indice++;
		}else{
			new_indice[i] = new_indice[representant[i]];
		}
	}
	#ifdef DEBUG
	for(int i=0 ; i<G.nb_sommets ; i++){
		cout<<"Sommet "<<i<<" d'indice "<<new_indice[i]<<"\n";
	}
	#endif
//	cout<<"Non colorie : "<<cpt_non_colorie<<" Solution courante : "<<solution_courante<<"\n";
	GC.nb_sommets = solution_courante+cpt_non_colorie;
//	cout<<solution_courante<<" taille="<<GC.nb_sommets<<"\n";
	GC.matrice_adjacence.resize(GC.nb_sommets);
	for(int i=0 ; i<GC.nb_sommets ; i++){
		GC.matrice_adjacence[i].resize(GC.nb_sommets);
	}
	GC.sommets_voisins.resize(GC.nb_sommets);
	for(int i=0 ; i<G.nb_sommets-1 ; i++){if(representant[i] == i){
		for(int j=i+1 ; j<G.nb_sommets ; j++){if(representant[j] == j){
			GC.matrice_adjacence[new_indice[i]][new_indice[j]] = matrice_adjacence[i][j];
			GC.matrice_adjacence[new_indice[j]][new_indice[i]] = matrice_adjacence[i][j];
			if(GC.matrice_adjacence[new_indice[i]][new_indice[j]] == 1){
				GC.sommets_voisins[new_indice[j]].push_back(new_indice[i]);
				GC.sommets_voisins[new_indice[i]].push_back(new_indice[j]);
			}
		}}
	}}

	free(representant_couleur);

	return ;
}

void creer_graphe_GC2(C_Graphe &G, C_Graphe &GC, int *coloration_courante, int solution_courante, int *representant,int *matrice_adjacence){
	int cpt_non_colorie;
	cpt_non_colorie = 0;
	//On construit un sommet par groupe de couleurs pour remplacer les sommets déja coloriés en élisant un représentant
	representant[G.nb_sommets+2-1] = 1;

	//Vecteur tampon matrice_adjacence
	for(int i=0 ; i<G.nb_sommets ; i++){
		if(coloration_courante[i] == -1){cpt_non_colorie++;}
		representant[i] = -1;
		for(int j=i ; j<G.nb_sommets ; j++){
			matrice_adjacence[i*G.nb_sommets+j] = 0;
			matrice_adjacence[j*G.nb_sommets+i] = 0;
		}
	}

	int *representant_couleur;
	representant_couleur = (int *)malloc(sizeof(int)*(solution_courante+1));
	for(int i=1 ; i<=solution_courante ; i++){
		representant_couleur[i] = -1;
	}
	for(int i=0 ; i<G.nb_sommets ; i++){
		#ifdef DEBUG
		cout<<i<<" "<<coloration_courante[i]<<" "<<representant_couleur[coloration_courante[i]]<<"\n";
		#endif
		if(coloration_courante[i] != -1){
			if(representant_couleur[coloration_courante[i]] == -1){
				representant_couleur[coloration_courante[i]] = i;
				representant[i] = i;
			}else{
				representant[i] = representant_couleur[coloration_courante[i]];
			}
		}else{representant[i] = i;}
	}
	for(int i=0 ; i<G.nb_sommets ; i++){
		if(representant[i] != i){
			for(int j=0 ; j<G.nb_sommets ; j++){
				matrice_adjacence[i*G.nb_sommets+j] = -1;
				matrice_adjacence[j*G.nb_sommets+i] = -1;
			}
		}
	}
	#ifdef DEBUG
	cout<<"Representants couleur\n";
	for(int i=1 ; i<=solution_courante ; i++){
		cout<<"Couleur "<<i<<" représentée par "<<representant_couleur[i]<<"\n";
	}
	cout<<"Representants\n";
	for(int i=0 ; i<G.nb_sommets ; i++){
		cout<<"Sommet "<<i<<" représenté par "<<representant[i]<<"\n";
	}
	#endif
	//On change les vecteurs d'adjacence des représentants
	for(int i=0 ; i<G.nb_sommets ; i++){
		for(int j=i ; j<G.nb_sommets ; j++){
			if(G.adj[i*G.nb_sommets+j]){
				if(!matrice_adjacence[representant[i]*G.nb_sommets+representant[j]]){GC.nb_aretes++;}
				matrice_adjacence[representant[i]*G.nb_sommets+representant[j]] = 1;
				matrice_adjacence[representant[j]*G.nb_sommets+representant[i]] = 1;
			}
		}
	}
	//On crée finalement la matrice d'adjacence de GC
	GC.nb_sommets = solution_courante+cpt_non_colorie;
	GC.matrice_adjacence.resize(GC.nb_sommets);
	GC.adj = (bool *)malloc(sizeof(bool)*GC.nb_sommets*GC.nb_sommets);
//	GC.sommets_voisins_bis = (int *)malloc(sizeof(int)*GC.nb_sommets*GC.nb_sommets);
//	GC.degre = (int *)malloc(sizeof(int)*GC.nb_sommets);
//	for(int i=0 ; i<GC.nb_sommets ; i++){GC.degre[i] = 0;}
	int cpt_i = 0;
	for(int i=0 ; i<G.nb_sommets-1 ; i++){if(representant[i] == i){
		int cpt_j = cpt_i+1;
		for(int j=i+1 ; j<G.nb_sommets ; j++){if(representant[j] == j){
			GC.adj[cpt_i*GC.nb_sommets+cpt_j] = (matrice_adjacence[i*G.nb_sommets+j]== 1?true:false);
			GC.adj[cpt_j*GC.nb_sommets+cpt_i] = (matrice_adjacence[i*G.nb_sommets+j]== 1?true:false);
//			if(GC.adj[cpt_i*GC.nb_sommets+cpt_j]){
//				GC.sommets_voisins_bis[cpt_j*GC.nb_sommets+GC.degre[cpt_j]] = cpt_i;
//				GC.degre[cpt_j]++;
//				GC.sommets_voisins_bis[cpt_i*GC.nb_sommets+GC.degre[cpt_i]] = cpt_j;
//				GC.degre[cpt_i]++;
//			}
			cpt_j++;
		}}
		cpt_i++;
	}}

	free(representant_couleur);
	return ;
}

void C_Graphe::maximum_cardinality_search(vector<int> &alpha, vector<int> &alpha_1){
	vector < int > set_size;//array of the size of each set int double array set
	vector < vector<int> > set;
	vector <int> size;

	set.resize(nb_sommets);
	set_size.resize(nb_sommets);
	size.resize(nb_sommets);
	alpha.resize(nb_sommets);
	alpha_1.resize(nb_sommets);
	
	for (int i=0 ; i<nb_sommets ; i++){
		set[i].resize(nb_sommets);
		set_size[i] = 0;
		for (int j=0 ; j<nb_sommets ; j++){
			set[i][j] = 0;
		}
	}
	int degre_max = 0, i_degre_max = -1;
	for (int i=0 ; i<nb_sommets ; i++){
		size[i] = 0;
		set[0][i] = 1;
		set_size[0]++;
		if (sommets_voisins[i].size() > degre_max){
			degre_max = sommets_voisins[i].size();
			i_degre_max = i;
		}
	}

	int i=nb_sommets-1, j=0, v;
	bool debut = true;
	v = i_degre_max;
	while(i>=0){
		//Choix du sommet a numeroter
		if (!debut){
			v=0;
			while(v<nb_sommets && set[j][v] == 0){v++;}
		}else{
			debut = false;
		}
		set[j][v] = 0;
		set_size[j]--;
		#ifdef DEBUG_chordal
		cout<<"Sommet choisi "<<v<<"\n";
		#endif

		//Numerotation
		alpha[v] = i; alpha_1[i] = v; size[v] = -1;

		//Update set et size
		for (list<int>::iterator w = sommets_voisins[v].begin() ; w != sommets_voisins[v].end() ; w++){
			if (size[*w] >= 0){
				set[size[*w]][*w] = 0;
				set_size[size[*w]]--;

				size[*w]++;

				set[size[*w]][*w] = 1;
				set_size[size[*w]]++;
				#ifdef DEBUG_chordal
				cout<<"Update size["<<*w<<"]="<<size[*w]<<"\t";
				cout<<"Update set["<<size[*w] - 1<<"]\t";
				cout<<"Update set["<<size[*w]<<"]\n";
				#endif
			}
		}
		i--;
		j++;
		while(j>=0 && set_size[j] == 0){j--;}
	}
	return;
}

bool C_Graphe::chordal(vector<int> &alpha, vector<int> &alpha_1){
	#ifdef DEBUG_chordal
	cout<<"Chordal test\n";
	#endif
	vector<int> index;
	vector<int> f;
	index.resize(nb_sommets);
	f.resize(nb_sommets);

	for(int i=0 ; i<nb_sommets ; i++){
		index[i] = -1;
		f[i] = -1;
	}

	#ifdef DEBUG_chordal
	for(int i=0 ; i<nb_sommets ; i++){
		cout<<"alpha["<<i<<"]="<<alpha[i]<<"\n";
	}
	for(int i=0 ; i<nb_sommets ; i++){
		cout<<"Sommet "<<i<<" de alpha "<<alpha_1[i]<<"\n";
	}
	#endif

	int w;
	for(int i=0 ; i<nb_sommets ; i++){
		w = alpha[i]; f[w] = w; index[w] = alpha[w];
		#ifdef DEBUG_chordal
		cout<<"Sommet "<<w<<"\n";
		#endif

		cout<<"Voisins\n";
		for(list<int>::iterator v = sommets_voisins[w].begin() ; v != sommets_voisins[w].end() ; v++){
			cout<<"Voisin "<<*v<<"\n";
			if (alpha[*v] < alpha[w]){
				index[*v] = alpha[w];
	 			if(f[*v] == *v){
					f[*v] = w;
				}
			}
		}

		#ifdef DEBUG_chordal
		cout<<"Followers\n";
		for(int i=0 ; i<nb_sommets ; i++){
			cout<<"f["<<i<<"]="<<f[i]<<"\n";
		}
		cout<<"Indexes\n";
		for(int i=0 ; i<nb_sommets ; i++){
			cout<<"index["<<i<<"]="<<index[i]<<"\n";
		}
		#endif 

		cout<<"Voisins\n";
		for(list<int>::iterator v = sommets_voisins[w].begin() ; v != sommets_voisins[w].end() ; v++){
			cout<<"Voisin "<<*v<<"\n";
			if (alpha[*v] < alpha[w] && index[f[*v]]<alpha[w]){
				cout<<"Probleme "<<*v<<" "<<f[*v]<<" "<<index[f[*v]]<<"\n";
				return false;
			}
		}
	}
	return true;
}

bool C_Graphe::comparability(){
	int k=0;
	FLAG=0;

	CLASS.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		CLASS[i].resize(nb_sommets);
		for(int j=0 ; j<nb_sommets ; j++){
			CLASS[i][j] = nb_sommets*2;//Pour class, undefined = nb_sommets*2
		}
	}

	for(int i=0 ; i<nb_sommets-1 ; i++){
		for(int j=i+1 ; j<nb_sommets ; j++){
			if (matrice_adjacence[i][j] == 1 && CLASS[i][j] == nb_sommets*2){
				k++;
				CLASS[i][j] = k;
				CLASS[j][i] = -k;
				explore(i,j,k);
				if (FLAG == 1){
					return false;
				}
			}
		}
	}
	return true;
}

void C_Graphe::explore(int i, int j, int k){
	if(FLAG == 1){return;}
	for (list<int>::iterator m = sommets_voisins[i].begin() ; m != sommets_voisins[i].end() ; m++){
		if (matrice_adjacence[j][*m] == 0 || ( CLASS[j][*m] != nb_sommets*2 && (abs(CLASS[j][*m])<k) ) ){
			if (CLASS[i][*m] == nb_sommets*2){
				CLASS[i][*m] = k;
				CLASS[*m][i] = -k;
				explore(i,*m,k);
			}else{if(CLASS[i][*m] == -k){
				FLAG = 1;
				return ;
			}}
		}
	}

	for (list<int>::iterator m = sommets_voisins[j].begin() ; m != sommets_voisins[j].end() ; m++){
		if (matrice_adjacence[i][*m] == 0 || ( CLASS[i][*m] != nb_sommets*2 && (abs(CLASS[i][*m])<k) ) ){
			if (CLASS[*m][j] == nb_sommets*2){
				CLASS[*m][j] = k;
				CLASS[j][*m] = -k;
				explore(*m,j,k);
			}else{if(CLASS[*m][j] == -k){
				FLAG = 1;
				return ;
			}}
		}
	}
	return ;
}

void C_Graphe::sommet_contraction(int u, int v){
	//Arbitrairement, on choisit de mettre tous les voisins de v dans u
	existe[v] = 0;
	nb_sommets--;
	for (list<int>::iterator voisin_v = sommets_voisins[v].begin() ; voisin_v != sommets_voisins[v].end() ; voisin_v++){
		if (existe[*voisin_v] == 1){
			if (matrice_adjacence[u][*voisin_v] == 0){
				nb_aretes++;
				matrice_adjacence[u][*voisin_v] = 1;
				matrice_adjacence[*voisin_v][u] = 1;
				sommets_voisins[u].push_back(*voisin_v);
				sommets_voisins[*voisin_v].push_back(u);
			}
		}
	}
	return;
}

void C_Graphe::aretes_ajout(int u, int v){
	matrice_adjacence[u][v] = 1;
	matrice_adjacence[v][u] = 1;
	sommets_voisins[u].push_back(v);
	sommets_voisins[v].push_back(u);
	nb_aretes++;
	return;
}

void C_Graphe::copier(C_Graphe &G){
	nb_sommets = G.nb_sommets;
	nb_aretes = G.nb_aretes;
	seed = G.seed;

	//Matrice d'adjacence et existe
	matrice_adjacence.resize(nb_sommets);
	existe.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		matrice_adjacence[i].resize(nb_sommets);
		existe[i] = G.existe[i];
		for (int j=0 ; j<nb_sommets ; j++){
			matrice_adjacence[i][j] = G.matrice_adjacence[i][j];
		}
	}

	//Liste de sommets
	sommets_voisins.resize(nb_sommets);
	for(int i=0 ; i<nb_sommets ; i++){
		for (list<int>::iterator voisin = sommets_voisins[i].begin() ; voisin != sommets_voisins[i].end() ; voisin++){
			sommets_voisins[i].push_back(*voisin);
		}
	}
	return;
}

C_Graphe::~C_Graphe(){
}
