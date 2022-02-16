
#include <ilcplex/ilocplex.h>

#include "dimacs.hpp"
#include "graph.hpp"
#include "options.hpp"

ILOSTLBEGIN

/************************************************************************************
 ************************************************************************************
 * Main
 ************************************************************************************
 ***********************************************************************************/

using weight = int;

int main(int argc, char* argv[])
{

    auto options = gc::parse(argc, argv);

    using dense_graph = gc::graph<gc::bitset>;
    dense_graph g;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = dense_graph{nv}; },
        [&](int u, int v) {
            if (u != v) {
                g.add_edge(u - 1, v - 1);
            }
        },
        [&](int, weight) {});

    auto sol{gc::brelaz_color(g)};
    int ub{*max_element(begin(sol), end(sol)) + 1};

    gc::clique_finder<gc::bitset> cf{g};
    cf.find_cliques(g.nodes, 0, g.size());

    IloEnv env;
    IloModel mod(env, "coloring");

    std::cout << "create variables\n";

    IloArray<IloIntVarArray> colors(env, g.nodes.size());
    IloIntVarArray usedColors(env, ub, 0, 1);
    IloExpr obj(env);

    std::cout << "used colors variables\n";
    for (int i = 0; i < ub; ++i) {
        mod.add(usedColors[i]);
        obj += usedColors[i];
        if (i < ub - 1)
            mod.add(usedColors[i] >= usedColors[i+1]);
    }

    std::cout << "node colors variables\n";
    for (auto v : g.nodes) {
        IloIntVarArray possible_colors(env, ub, 0, 1);
        colors[v] = possible_colors;
        for (int i = 0; i < ub; ++i) {
            mod.add(colors[v][i]);
        }
    }

    std::cout << "add constraints\n";
    for (auto u : g.nodes) {
        for (auto v : g.nodes) {
            if (u < v && g.matrix[u].contain(v)) {
                for (int i = 0; i < ub; ++i) {
                    IloExpr edge(env);
                    edge += colors[u][i];
                    edge += colors[v][i];
                    mod.add(edge <= 1);
                }
            }
        }
        IloExpr atleast(env);
        for (int i = 0; i < ub; ++i) {
            mod.add(colors[u][i] <= usedColors[i]);
            atleast += colors[u][i];
        }
        mod.add(atleast >= 1);
    }

    // for(int i=0; i<cf.num_cliques; ++i) {
    //      for(int c=0; c<ub; ++c) {
    //              IloExpr clique(env);
    //              for(auto v : cf.cliques[i]) {
    //                      clique += colors[v][c];
    //              }
    //              mod.add( clique <= 1 );
    //      }
    // }

    std::cout << "add clique constraints\n";

    auto maxidx{std::distance(begin(cf.clique_sz),
        std::max_element(
            begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques))};

    int c{0};
    for (auto v : cf.cliques[maxidx]) {
        mod.add(colors[v][c++] >= 1);
    }

    std::cout << "add objective\n";
    mod.add(IloMinimize(env, obj));

    IloCplex cplex(mod);
    cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::RandomSeed, options.seed);
		
		std::cout << "set seed " << options.seed << std::endl;

    if (cplex.solve()) {
        cplex.out() << "Objective value: " << cplex.getObjValue() << endl;
    } else {
        cplex.out() << "No solution available" << endl;
    }

    // std::cout << "end of program" << std::endl;
    return 0;
}
