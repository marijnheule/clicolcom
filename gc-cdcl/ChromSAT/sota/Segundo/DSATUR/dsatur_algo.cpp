#include "dsatur_algo.h"
#include "src/bitset.hpp"
//#define DEBUG

#include <iomanip>

#define NOEUDS 1000000

#define activate_time
#define activate_heur

//#define show_UB
//#define show_appel

// #define DEBUG_END
// #define DEBUG_PREP

/*************************************************************************************************************************/
/*************************************************************************************************************************/

int tiebreak = 3;

int *reorder;

int *F_V1;
int *F_V2;
bool *F;
int *candidats;
int *s;
bool *adj;

int size_candidats;

int profondeur_glob;

std::mt19937 random_generator;


/*************************************************************************************************************************/

int DSATUR_::DSATUR_h(C_Graphe& GC)
{
    int* coloration_courante;
    int* dsat_courant;
    int* candidats;

    coloration_courante = (int*)malloc(GC.nb_sommets * sizeof(int));
    dsat_courant = (int*)malloc(GC.nb_sommets * sizeof(int));
    candidats = (int*)malloc(GC.nb_sommets * sizeof(int));

    for (int i = 0; i < GC.nb_sommets; i++) {
        coloration_courante[i] = -1;
        dsat_courant[i] = 0;
    }
    int cpt = 0, cpt_boucle = 0;
    int choix_coul, choix_som, solution_courante = 0, size_candidats = 0;
    while (cpt_boucle < GC.nb_sommets) {
        // Choix du sommet
        choix_som = -1;
        int choix_dsat = -1;
        for (int j = 0; j < G.nb_sommets; j++) {
            int i = reorder[j];

            if (coloration_courante[i] == -1) {
                if (dsat_courant[i] > choix_dsat) {
                    choix_som = i;
                    choix_dsat = dsat_courant[i];

                    size_candidats = 0;
                    candidats[size_candidats] = i;
                    size_candidats++;

                } else {
                    if (dsat_courant[i]
                        == choix_dsat) { // En cas d'égalité on le rajoute
                        candidats[size_candidats] = i;
                        size_candidats++;
                    }
                }
            }
        }
        int choix_pass = -1;
        for (int i = 0; i < size_candidats; i++) {
            cpt = 0;
            for (int j = 0; j < size_candidats; j++) {
                if (i != j && adj[candidats[i] * G.nb_sommets + candidats[j]]) {
                    cpt += same(candidats[i], candidats[j]);
                }
            }
            if (cpt > choix_pass) {
                choix_som = candidats[i];
                choix_pass = cpt;
            }
        }
        // Choix de la couleur
        choix_coul = -1;
        bool ajout;
        for (int i = 1; i <= solution_courante + 1; i++) {
            ajout = true;
            choix_coul = i;
            for (int k = 0; k < GC.degre[choix_som]; k++) {
                if (coloration_courante
                        [GC.sommets_voisins_bis[choix_som * GC.nb_sommets + k]]
                    == choix_coul) {
                    ajout = false;
                    break;
                }
            }
            if (ajout) {
                break;
            }
        }
        coloration_courante[choix_som] = choix_coul;
        solution_courante = max(solution_courante, choix_coul);
        // Mise à jour dsat
        for (int k = 0; k < GC.degre[choix_som]; k++) {
            ajout = true;
            int tmp = GC.sommets_voisins_bis[choix_som * GC.nb_sommets + k];
            for (int l = 0; l < GC.degre[tmp]; l++) {
                if (coloration_courante
                        [GC.sommets_voisins_bis[tmp * GC.nb_sommets + l]]
                    == choix_coul) {
                    ajout = false;
                    break;
                }
            }
            if (ajout) {
                dsat_courant[tmp]++;
            }
        }
        cpt_boucle++;

    }

    if (solution_courante < UB) {
        store_solution(solution_courante, coloration_courante);
    }

    // for (int i = 0; i < init_n; ++i)
    //             meilleure_coloration[i] = coloration_courante[vertex_map[i]]
    //             - 1;

    free(coloration_courante);
    free(dsat_courant);
    free(candidats);

    return solution_courante;
}

void DSATUR_::store_solution(
    const int solution_courante, const int* coloration_courante)
{

    int maxc{0};
    // start by coloring the residual vertices
    for (int i = 0; i < G.nb_sommets; ++i) {
        // std::cout << i << " " << vmap[i] << " " << vertex_map[vmap[i]] <<
        // std::endl;

        auto c{coloration_courante[i] - 1};
        if (c > maxc)
            maxc = c;

        assert(vertex_map[vmap[i]] == i);
        meilleure_coloration[vmap[i]] = c;
    }

    assert(maxc == (solution_courante - 1));

    // now the dominated vertices
    assert(domdeleted.size() == dom_map.size());
    for (int i = domdeleted.size(); --i >= 0 ;) {

        assert(vertex_map[domdeleted[i]] == dom_map[i]);
        assert(meilleure_coloration[dom_map[i]] >= 0);

        // std::cout << "sol[" << domdeleted[i] << "] = cur[" << dom_map[i] <<
        // "]\n";

        meilleure_coloration[domdeleted[i]] = meilleure_coloration[dom_map[i]];
    }

    // if (lowdegdeleted.size() > 0) {
    //
    //     gc::bitset colors(0, solution_courante - 1, gc::bitset::full);
    //
    //     // now the low degree vertices
    //     for (int i = lowdegdeleted.size(); --i >= 0;) {
    //         auto v{lowdegdeleted[i]};
    //
    //         // assert(vertex_map[v] == -1);
    //
    //         for (int k = 0; k < G.degre[v]; k++) {
    //             auto c{
    //                 coloration_courante[G.sommets_voisins_bis[v * G.nb_sommets
    //                     + k]]
    //                 - 1};
    //             colors.fast_remove(c);
    //         }
    //
    //         // assert(!(colors.empty()));
    //
    //         meilleure_coloration[lowdegdeleted[i]] = colors.min();
    //     }
    // }
}

void DSATUR_::DSATUR_preprocessing()
{
#ifdef DEBUG_PREP
    cout << "PREPROCESSING\n";
#endif
    int* deleted_vertices = (int*)malloc(sizeof(int) * G.nb_sommets);
    int new_nb_sommets = G.nb_sommets;
    for (int i = 0; i < G.nb_sommets; i++) {
        deleted_vertices[i] = 0;
        vertex_map[i] = -1;
    }
#ifdef DEBUG_PREP
    cout << "RECHERCHE DOMINE\n";
#endif
    // On retire potentiellement les dominated vertices
    for (int i = 0; i < G.nb_sommets; i++) {
        for (int j = i; j < G.nb_sommets; j++) {
            if (i != j && G.degre[j] > G.degre[i]) {
                bool dominated = true;
                for (int k = 0; k < G.degre[i]; k++) {
                    bool found = false;
                    for (int l = 0; l < G.degre[j]; l++) {
                        if (G.sommets_voisins_bis[i * G.nb_sommets + k]
                            == G.sommets_voisins_bis[j * G.nb_sommets + l]) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        dominated = false;
                        break;
                    }
                }
                if (dominated) {
                    if (deleted_vertices[i] == 0) {
                        new_nb_sommets--;
                        vertex_map[i] = j;
                        domdeleted.push_back(i);
                        dom_map.push_back(j);

#ifdef DEBUG_PREP
                        cout << "Sommet " << i << " : ";
                        for (int k = 0; k < G.degre[i]; k++) {
                            cout << G.sommets_voisins_bis[i * G.nb_sommets + k]
                                 << " ";
                        }
                        cout << "\n";
                        cout << "Dominé par le sommet " << j << " : ";
                        for (int k = 0; k < G.degre[j]; k++) {
                            cout << G.sommets_voisins_bis[j * G.nb_sommets + k]
                                 << " ";
                        }
                        cout << "\n";
#endif
                    }
                    deleted_vertices[i] = 1;
                }
            }
        }
    }
#ifdef DEBUG_PREP
    cout << "RECHERCHE LB\n";
#endif
    // On retire les sommets tels que d(v) < LB
    for (int i = 0; i < G.nb_sommets; i++) {
        if (G.degre[i] < (LB - 1)) {
            if (deleted_vertices[i] == 0) {
                new_nb_sommets--;
            }
            deleted_vertices[i] = 2;
            lowdegdeleted.push_back(i);
        }
    }

    cout << "[trace] " << domdeleted.size() << " dominated vertices "
         << lowdegdeleted.size() << " low degree vertices\n";

    if (new_nb_sommets < G.nb_sommets) {
#ifdef DEBUG_PREP
        cout << G.nb_sommets - new_nb_sommets << " sommet(s) a enlever\n";
#endif
        // On crée un nouveau graphe en enlevant tous ces sommets
        C_Graphe Gnew;
        Gnew.nb_sommets = new_nb_sommets;
        Gnew.nb_aretes = 0;

        Gnew.adj
            = (bool*)malloc(Gnew.nb_sommets * Gnew.nb_sommets * sizeof(bool));
        Gnew.sommets_voisins_bis
            = (int*)malloc(sizeof(int) * Gnew.nb_sommets * Gnew.nb_sommets);
        Gnew.degre = (int*)malloc(sizeof(int) * Gnew.nb_sommets);

        for (int i = 0; i < Gnew.nb_sommets; i++) {
            Gnew.degre[i] = 0;
            for (int j = 0; j < Gnew.nb_sommets; j++) {
                Gnew.sommets_voisins_bis[i * Gnew.nb_sommets + j] = 0;
                Gnew.adj[i * Gnew.nb_sommets + j] = false;
            }
        }

        int indice_i = 0;
        for (int i = 0; i < G.nb_sommets; i++) {
            if (deleted_vertices[i] == 0) {

                assert(vertex_map[i] == -1);
                vertex_map[i] = indice_i;
                assert(vmap.size() == indice_i);
                vmap.push_back(i);

                if (i < G.nb_sommets - 1) {
                    int indice_j = indice_i + 1;
                    for (int j = i + 1; j < G.nb_sommets; j++) {
                        if (deleted_vertices[j] == 0) {
                            Gnew.adj[indice_i * Gnew.nb_sommets + indice_j]
                                = G.adj[i * G.nb_sommets + j];
                            Gnew.adj[indice_j * Gnew.nb_sommets + indice_i]
                                = G.adj[i * G.nb_sommets + j];
                            if (G.adj[i * G.nb_sommets + j]) {
                                Gnew.nb_aretes++;

                                Gnew.sommets_voisins_bis[(indice_i
                                                             * Gnew.nb_sommets)
                                    + Gnew.degre[indice_i]]
                                    = indice_j;
                                Gnew.degre[indice_i]++;

                                Gnew.sommets_voisins_bis[(indice_j
                                                             * Gnew.nb_sommets)
                                    + Gnew.degre[indice_j]]
                                    = indice_i;
                                Gnew.degre[indice_j]++;
                            }
                            indice_j++;
                        }
                    }
                    indice_i++;
                }
            }
        }
        free(G.sommets_voisins_bis);
        free(G.degre);
        free(G.adj);
        G = Gnew;
    } else {
        for (int i = 0; i < G.nb_sommets; i++) {
            vertex_map[i] = i;
            vmap.push_back(i);
        }
    }
    free(deleted_vertices);
    return;
}

int DSATUR_::DSATUR_algo(C_Graphe& G_param, double time_param, int regle_input,
    int LB_input, int UB_input, int seed)
{
    // #ifdef DEBUG
    cout << "[options] Paramètres :\n";
		cout << "[options] seed=" << seed << "\n";
    cout << "[options] time limit=" << time_param << "\n";
    // cout<<"[options] borne="<<borne_input<<"\n";
    // cout<<"[options] regle="<<regle_input<<"\n";
    // cout<<"[options]
    // profondeur=["<<borne_inf_prof_input<<","<<borne_sup_prof_input<<"]\n";
    // cout<<"[options]
    // gap=["<<borne_inf_gap_input<<","<<borne_sup_gap_input<<"]\n";
    // #endif

    cout << "[statistics] lb ub time conflicts\n";

		random_generator.seed(seed);

    start = clock();

    // Initialisation parametres
    G = G_param;
    LB = max(LB_input, 2);
    UB = min(UB_input, G.nb_sommets);

    cout << "[trace] preprocessing at "
         << (double)(clock() - start) / CLOCKS_PER_SEC << "\n";

    // Initialisation
    init_n = G.nb_sommets;
    vertex_map = (int*)malloc(sizeof(int) * init_n);

    DSATUR_preprocessing();

    time_limit = time_param;
    regle_sel = regle_input;

    quit = -1;
    nombre_noeuds = 1;

    cout << "[trace] alloc at " << (double)(clock() - start) / CLOCKS_PER_SEC
         << "\n";

    meilleure_coloration = (int*)malloc(sizeof(int) * init_n);
    coloration_courante = (int*)malloc(sizeof(int) * G.nb_sommets);
    dsat_courant = (int*)malloc(sizeof(int) * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets; i++) {
        meilleure_coloration[i] = -1;
        coloration_courante[i] = -1;
        dsat_courant[i] = 0;
    }
    solution_courante = 0;

    reorder = (int*)malloc(sizeof(int) * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets; i++) {
        reorder[i] = i;
    }
    candidats = (int*)malloc(sizeof(int) * G.nb_sommets);
    F_V1 = (int*)malloc(sizeof(int) * G.nb_sommets);
    F_V2 = (int*)malloc(sizeof(int) * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets; i++) {
        F_V1[i] = 1;
        F_V2[i] = 1;
    }
    F = (bool*)malloc(sizeof(bool) * G.nb_sommets * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets * G.nb_sommets; i++) {
        F[i] = true;
    }
    s = (int*)malloc(sizeof(int) * G.nb_sommets * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets * G.nb_sommets; i++) {
        s[i] = 1;
    }
    size_candidats = G.nb_sommets;
    adj = (bool*)malloc(sizeof(bool) * G.nb_sommets * G.nb_sommets);
    for (int i = 0; i < G.nb_sommets; i++) {
        for (int j = 0; j < G.nb_sommets; j++) {
            adj[i * G.nb_sommets + j] = G.adj[i * G.nb_sommets + j];
        }
    }

    cout << "[trace] upper bound at "
         << (double)(clock() - start) / CLOCKS_PER_SEC << "\n";

    // #ifdef activate_heur
    int pUB = DSATUR_h(G);
    if (pUB < UB)
        UB = pUB;

    if (reduction != NULL) {
        pUB = reduction->extend_solution(meilleure_coloration, UB);
    }

    std::cout << "[data] lb = " << std::setw(4) << std::left << LB
              << "| ub = " << std::setw(4) << std::left << pUB
              << "| time = " << std::setw(10) << std::left
              << std::setprecision(4)
              << (double)(clock() - start) / CLOCKS_PER_SEC
              << "| conflicts = " << std::setw(10) << std::left << nombre_noeuds
              << std::endl;

    // #endif

    cout << "[trace] lower bound at "
         << (double)(clock() - start) / CLOCKS_PER_SEC << "\n";

    // Calcul de la borne inf
    int* solution = (int*)malloc(sizeof(int) * G.nb_sommets);
    int size_solution = 0;
    LB = clique_cliquer_init(G, 0, solution, &size_solution);

    std::cout << "[data] lb = " << std::setw(4) << std::left << LB
              << "| ub = " << std::setw(4) << std::left << pUB
              << "| time = " << std::setw(10) << std::left
              << std::setprecision(4)
              << (double)(clock() - start) / CLOCKS_PER_SEC
              << "| conflicts = " << std::setw(10) << std::left << nombre_noeuds
              << std::endl;

    for (int i = 0; i < size_solution; i++) {
        int tmp = reorder[i];
        reorder[i] = reorder[solution[i]];
        reorder[solution[i]] = tmp;
    }
    for (int i = 0; i < size_solution; i++) {
        coloration_courante[reorder[i]] = i + 1;
        mise_a_jour_dsat_noupdate(reorder[i], i + 1);
        solution_courante++;
    }
    for (int i = size_solution; i < G.nb_sommets; i++) {
        int max = -1;
        int i_max = -1;
        for (int j = i; j < G.nb_sommets; j++) {
            if (G.degre[reorder[j]] > max) {
                max = G.degre[reorder[j]];
                i_max = j;
            }
        }
        int tmp = reorder[i];
        reorder[i] = reorder[i_max];
        reorder[i_max] = tmp;
    }
    free(solution);

#ifdef DEBUG
	cout<<"LB="<<LB<<" "<<"UB="<<UB<<"\n";
#endif
	
	cout << "[trace] search at " << (double)(clock() - start)/CLOCKS_PER_SEC << "\n";
	
	
	profondeur_glob = 0;
        if (DSATUR_algo_rec(solution_courante))
            LB = UB;
        time_spent = (double)(clock() - start) / CLOCKS_PER_SEC;

        // Free variables globales
        free(reorder);
        free(candidats);

        free(G.sommets_voisins_bis);
        free(G.degre);
        free(G.adj);

        free(F);
        free(s);
        free(adj);

        return UB;
}

bool DSATUR_::DSATUR_algo_rec(int profondeur)
{
    profondeur_glob = profondeur;

#ifdef DEBUG_END
    cout << "profondeur " << profondeur << ": " << solution_courante << " ["
         << LB << ".." << UB << "]"
         << "\n";
#endif
#ifdef activate_time
	double time = (double)(clock() - start)/CLOCKS_PER_SEC;
#endif
	if(UB <= LB){//Solution optimale trouvée, on arrête tout
#ifdef DEBUG_END
            cout << "OPT PROUVE\n";
#endif
            quit = 1;
            return true;
        }
#ifdef activate_time
	if(time >= time_limit){//Timer
#ifdef DEBUG_END
            cout << "TIMEOUT\n";
#endif
            quit = 1;
            return false;
        }
#endif
	if(quit == 1){//Arret récursion car la solution optimale a été trouvée
#ifdef DEBUG_END
            cout << "QUIT\n";
#endif
            return true;
        }
        if (profondeur == G.nb_sommets) { // Feuille
            nombre_noeuds++;

#ifdef DEBUG
            cout << "FEUILLE\n";
            afficher_solution_courante(profondeur);
#endif
            if (solution_courante < UB) {
#ifdef DEBUG
                cout << "Amélioration solution courante : " << solution_courante
                     << "\n";
#endif
#ifdef show_UB
#ifdef activate_time
                cout << "Amélioration solution courante : " << solution_courante
                     << " " << nombre_noeuds << " " << time << "\n";
#else
                cout << "Amélioration solution courante : " << solution_courante
                     << " " << nombre_noeuds << "\n";
#endif
#endif
                UB = solution_courante;
                store_solution(solution_courante, coloration_courante);

                // if(print_progress)
                int pUB = UB;
                if (reduction != NULL) {
                    pUB = reduction->extend_solution(meilleure_coloration, UB);
                }
                std::cout << "[data] lb = " << std::setw(4) << std::left << LB
                          << "| ub = " << std::setw(4) << std::left << pUB
                          << "| time = " << std::setw(10) << std::left
                          << std::setprecision(4)
                          << (double)(clock() - start) / CLOCKS_PER_SEC
                          << "| conflicts = " << std::setw(10) << std::left
                          << nombre_noeuds << std::endl;
            }
//			}
		#ifdef DEBUG
		cout<<"FIN FEUILLE\n";
		#endif
                return true;
    } else { // Pas une feuille
        nombre_noeuds++;
#ifdef DEBUG
        cout << "ELSE\n";
        afficher_solution_courante(profondeur);
#endif

        // Selection sommet
        int choix_sommet;
        int mu;
        choix_sommet = choisir_sommet(&size_candidats); // if(regle_sel == 1)
        if (regle_sel == 2) {
            mu = solution_courante;
            for (int i = 0; i < size_candidats; i++) {
#ifdef DEBUG
                cout << "Candidat : " << candidats[i] << "\n";
#endif
                mu = mu - dsat_courant[candidats[i]];
            }
            if (mu <= tiebreak) {
                choix_sommet = choisir_sommet_PASS(size_candidats);
            }
        }
#ifdef DEBUG
        cout << "Sommet choisi : " << choix_sommet << "\n";
#endif

        // Save des anciens paramètres/changements
        int old_couleur = coloration_courante[choix_sommet];
        int old_solution_courante = solution_courante;

        int* liste_changements_F;
        int size_changements_F;

        liste_changements_F = (int*)malloc(sizeof(int) * G.nb_sommets);

        // Recursion
        for (int i = 1; i <= solution_courante + 1; i++) {
            if (F[choix_sommet * G.nb_sommets + i]) {
#ifdef DEBUG
                cout << "Couleur choisie : " << i << "\n";
#endif
                // Si la couleur est nouvelle, on le signale

                // Mise a jour solution courante
                solution_courante = max(solution_courante, i);

                if (solution_courante < UB) {
                    // Mise a jour DSAT et F
                    mise_a_jour_dsat(choix_sommet, i, liste_changements_F,
                        &size_changements_F);

                    // Mise a jour coloration
                    coloration_courante[choix_sommet] = i;

                    // Appel recursif
                    DSATUR_algo_rec(profondeur + 1);

                    // Demise a jour coloration
                    coloration_courante[choix_sommet] = old_couleur;

                    // Demise a jour solution courante
                    solution_courante = old_solution_courante;

                    // Demise a jour DSAT et F
                    F[choix_sommet * G.nb_sommets + i] = true;
                    for (int j = 0; j < size_changements_F; j++) {
                        F[liste_changements_F[j] * G.nb_sommets + i] = true;
                        dsat_courant[liste_changements_F[j]]--;
                    }
                } else {
                    solution_courante = old_solution_courante;
                }
            }
        }
#ifdef DEBUG
        cout << "FIN ELSE\n";
#endif
        free(liste_changements_F);
    }
    return true;
}

int DSATUR_::choisir_sommet(int* size_candidats)
{
    int choix_sommet = -1;
    int choix_dsat = -1;
    (*size_candidats) = 0;
    for (int j = 0; j < G.nb_sommets; j++) {
        int i = reorder[j];

        if (coloration_courante[i] == -1) {
            if (dsat_courant[i] > choix_dsat) {
                choix_sommet = i;
                choix_dsat = dsat_courant[i];

                (*size_candidats) = 0;
                candidats[(*size_candidats)] = i;
                (*size_candidats)++;

            } else {
                if (dsat_courant[i]
                    == choix_dsat) { // En cas d'égalité on le rajoute
                    candidats[(*size_candidats)] = i;
                    (*size_candidats)++;
                }
            }
        }
    }

    return choix_sommet;
}

int DSATUR_::choisir_sommet_PASS(int size_candidats){
	int choix_sommet = -1;
	int choix_pass = -1;
	int cpt;
	int nb_best = 0;
	for(int i=0 ; i<size_candidats ; i++){
		cpt = 0;
		for(int j=0 ; j<size_candidats ; j++){if(i != j && adj[candidats[i]*G.nb_sommets+candidats[j]]){
			cpt += same(candidats[i],candidats[j]);
		}}
		if(cpt > choix_pass){
			choix_sommet = candidats[i];
			choix_pass = cpt;
			nb_best = 1;
		} 
		else if(cpt == choix_pass){
			// cout << "random";
			++nb_best;
			if(random_generator() % nb_best == 0) {
				choix_sommet = candidats[i];
				// cout << "*";
			}
			// cout << endl;
		}
	}
	return choix_sommet;
}

int DSATUR_::same(int v1, int v2){
	#ifdef DEBUG
	cout<<"SAME "<<v1<<" "<<v2<<"\n";
	#endif
	int cpt = 0;

	for(int i=0 ; i<solution_courante+1 ; i++){
		if(F[v1*G.nb_sommets+i] && F[v2*G.nb_sommets+i]){
			cpt++;
		}
	}

	return cpt;
}

void DSATUR_::afficher_solution_courante(int profondeur){
	cout<<"Coloration courante\n";
	cout<<"LB="<<LB<<" UB="<<UB<<" profondeur="<<profondeur<<"\n";
	cout<<"Coloration :\n";
	for(int i=0 ; i<G.nb_sommets ; i++){
		cout<<"Sommet "<<i<<" de couleur "<<coloration_courante[i]<<" et de DSAT "<<dsat_courant[i]<<"\n";
	}
	return;
}

void DSATUR_::choisir_couleurs(int choix_sommet, int *choix_couleurs, int *size_couleurs){
	#ifdef DEBUG
	cout<<"Debut choisir couleur\n";
	#endif
	(*size_couleurs) = 0;

	for(int i=1 ; i<=solution_courante ; i++){
		#ifdef DEBUG
		cout<<"Couleur candidate : "<<i<<"\n";
		#endif
		if(F[choix_sommet*G.nb_sommets+i]){//Si la couleur est disponible, alosr on la rajoute
			choix_couleurs[(*size_couleurs)] = i;
			(*size_couleurs)++;
		}
	}
	choix_couleurs[(*size_couleurs)] = solution_courante+1;
	(*size_couleurs)++;
	#ifdef DEBUG
	cout<<"Fin choisir couleur\n";
	#endif
	return ;
}

void DSATUR_::mise_a_jour_dsat(int choix_sommet, int choix_couleur, int *liste_changements_F, int *size_changements_F){
	(*size_changements_F) = 0;

	F[choix_sommet*G.nb_sommets+choix_couleur] = false;//La couleur qu'on vient d'assigner n'est par définition plus disponible
	
	for(int i=0 ; i<G.degre[choix_sommet] ; i++){
		if(F[G.sommets_voisins_bis[choix_sommet*G.nb_sommets+i]*G.nb_sommets+choix_couleur]){//Si la couleur était disponible, on la rend indisponible
			int tmp = G.sommets_voisins_bis[choix_sommet*G.nb_sommets+i];
			//Changement F
			F[tmp*G.nb_sommets+choix_couleur] = false;
			liste_changements_F[(*size_changements_F)] = tmp;
			(*size_changements_F)++;
			//Changement dsat
			dsat_courant[tmp]++;
		}
	}
	return;
}

void DSATUR_::mise_a_jour_dsat_noupdate(int choix_sommet, int choix_couleur){

	F[choix_sommet*G.nb_sommets+choix_couleur] = false;//La couleur qu'on vient d'assigner n'est par définition plus disponible
	
	for(int i=0 ; i<G.degre[choix_sommet] ; i++){
		if(F[G.sommets_voisins_bis[choix_sommet*G.nb_sommets+i]*G.nb_sommets+choix_couleur]){//Si la couleur était disponible, on la rend indisponible
		//Changement F
			F[G.sommets_voisins_bis[choix_sommet*G.nb_sommets+i]*G.nb_sommets+choix_couleur] = false;
		//Changement dsat
			dsat_courant[G.sommets_voisins_bis[choix_sommet*G.nb_sommets+i]]++;
		}
	}
	return;
}

DSATUR_::~DSATUR_(){
    free(meilleure_coloration);
    free(coloration_courante);
    free(dsat_courant);
}
