#ifndef __CG_SPARSE_GRAPH_HH
#define __CG_SPARSE_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"
#include "vertices_vec.hpp"

#include <algorithm>
#include <list>
#include <vector>

namespace gc
{

/*****************************************************************************/
/*					CLASS FOR SPARSE GRAPH

using adjacency_OLD vector describing neighbors of each vertex				 */
/*****************************************************************************/

class sparse_graph
{
public:
	IntStack nodes;
	vertices_vec nodeset;

	std::vector<std::vector<int>> adjacency_OLD; // To be given up
	
	std::vector<vertices_vec> adjacency;
	
	// Constructors
    sparse_graph() {}
	//~sparse_graph() {clear();}
	
	// With nv vertices, empty but reserved adjacency_OLD
    explicit sparse_graph(int nv)
    {
        nodes.reserve(nv);
        nodes.fill(); // Set nodes.size_ at nodes.list_.size()
		adjacency.reserve(nv); // Vector of size nv of empty vertices_vec, no neighboors yet
		vertices_vec temp;
		for (int i = 0; i < nv; ++i){
			adjacency.push_back(temp);
		}
    }
	
    sparse_graph(sparse_graph&) = default;
    sparse_graph(sparse_graph&&) = default;
    sparse_graph& operator=(const sparse_graph&) = default;
    sparse_graph& operator=(sparse_graph&&) = default;
	
	int capacity() const { return adjacency.size(); }
	
	int num_edges() const
	{
		int num_e(0);
		for(std::vector<vertices_vec>::const_iterator it = adjacency.begin(); it != adjacency.end(); ++it)
			num_e += (*it).size();
		num_e /= 2;
		return num_e;
	}

	int sparsity() const { return 2*num_edges()/(capacity()*(capacity()-1)); }
	
	// To be used in the case u > v (Due to the ordering of edges description in the graph files (DIMACS)
	void safe_add_edge(const int u, const int v)
	{	 
		std::vector<int>::iterator it;
		it = find(adjacency[v].vertices.begin(), adjacency[v].vertices.end(), u);
  		if (it == adjacency[v].vertices.end()) add_edge(u,v);// u is not already in adjacency[v]	 
	}

	void add_edge(int u, int v)
	{
		adjacency[u].push_back(v);
		adjacency[v].push_back(u);
	}

	// Safe ? May not be efficient 
	void remove_node(int u)
	{
		for (auto v : adjacency[u]) remove_edge(u, v);
		adjacency.erase(adjacency.begin() + u);	
	}

	// May be optimized (find + erase in vector -> moving left elements)
	// WARING: Only removing u from neighbors of v as u is to be removed from the graph anyway
	void remove_edge(int u, int v)
	{
		/*
		auto found = std::find(adjacency[v].begin(); adjacency[v].end(), u)
		if (found != adjacency[v].end()) adjacency[v].erase(found); 
		*/
		// else report error
		adjacency[v].remove(u);
	}

    void clear()
    {
        for (auto v : nodes) {
            adjacency[v].vertices.clear();
        }
		adjacency.clear();
        nodes.clear();
    }

	void describe(std::ostream& os) const;

	std::ostream& display_adjacency(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const sparse_graph& x);


/****************************************************************/
//					SPARSE CLIQUE FINDER						//
/****************************************************************/

struct sparse_clique_finder {
	const sparse_graph& g;
    std::vector<vertices_vec> cliques;
    std::vector<int> clique_sz;
    std::vector<vertices_vec> candidates; // Neighbors of the clique members
    std::vector<int> last_clique;
    int num_cliques;

    sparse_clique_finder(const sparse_graph& g);
	
		// clear previously cached results
		void clear();
		// initialize a new clique
		void new_clique();
    // initialize a new color
    void new_color();
    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);
    // insert v into the col^th color. assumes it fits. Puts vertices
    // added from candidates[i] into diff
    // void insert_color(int v, int col, bitset& diff);
    // heuristically find a set of cliques and return the size of the
    // largest
		

    template <class ordering> int find_cliques(ordering o, const int limit=0xfffffff) // BUGGY
    {
        clear();
        if (o.size() == 0)
            return 0;
				
		int max = 1;
				
        for (auto u : o) {
				
			//if(g.adjacency[u].size() < max) continue; // Put aside nodes with degree < size of the biggest clique so far
					
            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (std::find(candidates[i].vertices.begin(), candidates[i].vertices.end(), u) != candidates[i].vertices.end()) {
                    found = true;
                    insert(u, i); std::cout << "inserting vertex " << u+1 << " in clique " << i+1 << std::endl;
						
					if(cliques[i].size() > max) {
					max = cliques[i].size();
					}
                }

            if (!found && num_cliques < limit) {
                new_clique(); std::cout << "New clique " << num_cliques << std::endl;
                insert(u, num_cliques - 1); std::cout << "inserting vertex " << u+1 << " in NEW clique " << num_cliques << std::endl;
            }
        }

		// Second loop , could iterate over left candidates of cliques instead of all the vertices ? Too expensive.
        for (auto u : o) {
            for (int i = last_clique[u] + 1; i < num_cliques; ++i)
                if (std::find(candidates[i].vertices.begin(), candidates[i].vertices.end(), u) != candidates[i].vertices.end()) {
                    insert(u, i); std::cout << "inserting vertex " << u+1 << " in clique " << i+1 << std::endl;
                }
        }

        return *std::max_element(
            begin(clique_sz), begin(clique_sz) + num_cliques);
    }
	
	void display();
};


/****************************************************************/
//					Minimum degree elimination game				//
/****************************************************************/
struct Min_deg_elim_game {
	const sparse_graph& g;
	sparse_graph g_filled;
	sparse_graph g_section;

	std::vector<int> reverse_ordering; 

	Min_deg_elim_game(const sparse_graph& g);
	void elimination_game(std::vector<int> ordering);
		
};

/****************************************************************/
//					BRONKERBOSCH CLIQUE FINDER					//
/****************************************************************/

// Finding all maximal cliques (recursive)
struct BronKerbosch {
	const sparse_graph& g;

	// Collected results
	std::vector<std::vector<int>> maximal_cliques;
	std::vector<int> clique_sz;
	int num_cliques;
	
	// Algo variables
	std::vector<int> actual_clique;
	std::vector<int> candidates;
	std::vector<int> banned;

	// Vertices ordering
	// Degree
	std::vector<int> by_degree;
	std::vector<int> degree;
	void order_by_degree();

	// Degeneracy
	

	// Constructor
	BronKerbosch(const sparse_graph& g);
	
	void clear();
	
	void add_max_clique(std::vector<int>& clique);
	
	void find_cliques_withoutPivot(std::vector<int> clique, std::vector<int> candidates, std::vector<int> banned);

	// TODO 
	void find_cliques_withPivot(std::vector<int> clique, std::vector<int> candidates, std::vector<int> banned);
	int find_pivot(std::vector<int> const& v1, std::vector<int> const& v2);

	void display();

	void bronkerbosch_calls_display(std::vector<int> const& clique, std::vector<int> const& candidates, std::vector<int> const& banned);

	// Intersection, union
	std::vector<int> eliminate(std::vector<int> const& v1,std::vector<int> const& v2);
	std::vector<int> intersect(std::vector<int> const & v1, std::vector<int> const & v2);
	void unite_vector_element_ref(std::vector<int> & v,const int i); // modifies v
	std::vector<int> unite_vector_element_val(std::vector<int> v,const int i); // return a copy of v with i inserted
};

} // namespace gc

#endif

