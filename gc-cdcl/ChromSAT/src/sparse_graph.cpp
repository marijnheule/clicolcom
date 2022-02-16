
#include "sparse_graph.hpp"
#include "utils.hpp"
#include <minicsp/mtl/Heap.h>
#include <algorithm>
#include <iomanip>

namespace gc
{

// DISPLAY
void sparse_graph::describe(std::ostream& os) const
{
    os << "# vertices = " << capacity() << std::endl;
}

std::ostream& sparse_graph::display_adjacency(std::ostream& os) const
{
	os << "Displaying adjacency\n";
	int i(1);
    for (std::vector<vertices_vec>::const_iterator it = adjacency.begin(); it < adjacency.end(); ++it) {
        os << "Vertex " << i << "   Neighbors : ( ";
		i++;
		for (std::vector<int>::const_iterator itN = (*it).vertices.begin(); itN < (*it).vertices.end(); ++itN) {
			os << *itN+1 << " ";
		}
	    os << ")";
		os << " degree =" << (*it).vertices.size() << "\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const sparse_graph& x)
{
    return x.display_adjacency(os);
}

void sparse_clique_finder::display()
{
	std::cout << " Maximal cliques " << std::endl;
	int i(1);
	for(std::vector<vertices_vec>::const_iterator it = cliques.begin(); it != cliques.end(); ++it) {
		if (!(*it).empty()){
			std::cout << " Clique " << i << " { ";
			i++;
			for (std::vector<int>::const_iterator itN = (*it).vertices.begin(); itN < (*it).vertices.end(); ++itN) {
				std::cout << *itN+1 << " ";
			}
			std::cout << " } size " << (*it).vertices.size() <<std::endl;
		}
	}
}

/****************************************************************/
//					SPARSE CLIQUE FINDER						//
/****************************************************************/
sparse_clique_finder::sparse_clique_finder(const sparse_graph& g)
    : g(g)
    , num_cliques(1)
{
	last_clique.resize(g.capacity()); // OK
    cliques.resize(g.capacity()); // Could we divide it by 2 ?
    clique_sz.resize(g.capacity());
    candidates.resize(g.capacity());
	// Initialization of the candidates of the first clique (all vertices from g.nodes)
	candidates[0].vertices.reserve(g.capacity()); // Allocate the necessary space
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates[0].push_back(*it);
	}
}

void sparse_clique_finder::clear() { num_cliques = 0; }

// Declare a new clique (it should be associated with a vertex to put in)
void sparse_clique_finder::new_clique()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].vertices.clear();
	// Why set its size to 0 ?
    clique_sz[num_cliques] = 0;
	// All vertices candidates ?
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates[num_cliques].push_back(*it);
	}    
    ++num_cliques;
}

void sparse_clique_finder::new_color()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].vertices.clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].vertices.clear();
    ++num_cliques;
}

void sparse_clique_finder::insert(int v, int clq)
{
    cliques[clq].push_back(v);
    ++clique_sz[clq];
	// HERE we want candidates[clq] to be the intersected with N(v)
	// Should be sorted first (OK by construction)
	//candidates[clq].safe_intersect_with(g.adjacency[v]); // const& sparse_graph g in sparse_clique_finder !!! Anyway no need to sort here
	candidates[clq].intersect_with(g.adjacency[v]); 

/*
	std::vector<int>::iterator it;
	it = set_intersection(candidates[clq].vertices.begin(), candidates[clq].vertices.end(), g.adjacency[v].vertices.begin(), g.adjacency[v].vertices.end() , candidates[clq].vertices.begin());
	candidates[clq].vertices.resize(it-candidates[clq].vertices.begin());
*/	
    last_clique[v] = clq;
}


/*
void sparse_clique_finder::insert_color(int v, int clq, bitset& diff)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    diff.copy(g.matrix[v]);
    diff.setminus_with(candidates[clq]);
    candidates[clq].union_with(g.matrix[v]);
}
*/

/****************************************************************/
//					Minimum degree elimination game				//
/****************************************************************/

Min_deg_elim_game::Min_deg_elim_game(const sparse_graph& g)
	:g(g)
{
	g_filled = g;
	g_section = g;

	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		reverse_ordering.push_back(*it);
	} 
}

void Min_deg_elim_game::elimination_game(std::vector<int> ordering)
{
	int num_fill(0);
	for (auto o : ordering) { // o is int
		// Check every pair n1n2 of neighbors of o
		for (auto n1 : g_filled.adjacency[o]) { // n1 si vertices_vec
			for(std::vector<int>::iterator n2 = g_filled.adjacency[o].begin() + 1; n2 != g_filled.adjacency[o].end(); ++n2) {
				// If n1n2 not neighors add to g_section and g_filled
				auto found = std::find(g_filled.adjacency[n1].begin(), g_filled.adjacency[n1].end(), *n2);
				if (found != std::end(g_filled.adjacency[n1])) {
					g_filled.add_edge(n1, *n2);
					g_section.add_edge(n1, *n2);
					++num_fill;
				}
				
			}
		}
		g_section.remove_node(o);
	}	
}

/****************************************************************/
//						BRONKERBOSCH							//
/****************************************************************/

BronKerbosch::BronKerbosch(const sparse_graph& g)
	:g(g)
	, num_cliques(0)
{
	maximal_cliques.resize(g.capacity()); // Look for an upper bound to #maximal_cliques
										  // to prevent memory issues
    clique_sz.reserve(g.capacity());
    //actual_clique.resize(g.capacity());
	candidates.reserve(g.capacity()); // 
	banned.reserve(g.capacity());	// Too large

	// Set the candidates as all nodes
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates.push_back(*it);
	}
}

void BronKerbosch::clear() { num_cliques = 0; }

// Add a clique to maximal_cliques
void BronKerbosch::add_max_clique(std::vector<int>& clique)
{
	std::cout << "								ADDING CLIQUE " << std::endl;
	//assert(num_cliques < g.capacity());
	if (num_cliques >= g.capacity())	
		maximal_cliques.resize(maximal_cliques.size()+1);
    maximal_cliques[num_cliques].clear();
	for(std::vector<int>::const_iterator it = clique.begin(); it < clique.end(); ++it ) {
		maximal_cliques[num_cliques].push_back(*it);
	}
	clique_sz[num_cliques] = maximal_cliques[num_cliques].size();    
    ++num_cliques;
}

std::vector<int> BronKerbosch::intersect(std::vector<int> const& v1,std::vector<int> const& v2) // v1 & v2 are assumed to be sorted
{
	std::vector<int> v_intersection;
	std::set_intersection(v1.begin(), v1.end(),
						  v2.begin(), v2.end(),
						  std::back_inserter(v_intersection));	
	return v_intersection;
}

// Return v1 \ v2 (Vector?)
// BUGGED
	// Assume v1 and v2 are sorted (?)
std::vector<int> BronKerbosch::eliminate(std::vector<int> const& v1,std::vector<int> const& v2)
{
	std::vector<int> v_elim;
	v_elim.reserve(v1.size());
	// Pour chaque element de v1, le mettre dans v_elim si il n'est pas dans v2
	for(std::vector<int>::const_iterator it1 = v1.begin(); it1 != v1.end(); ++it1 ){
		//std::cout << "Ok in for loop in eliminate" << std::endl;
		auto it2 = find(v2.begin(), v2.end(), *it1);
		std::cout << "*it1 = " << *it1 << "  *it2 = " << *it2 << std::endl;
  		if (it2 == v2.end()){
			std::cout << "Pushing in v_elim : " << *it1+1 << std::endl;
			v_elim.push_back(*it1);}}
	return v_elim;
}

// std::vector::insert ... (Vector?)
std::vector<int> BronKerbosch::unite_vector_element_val(std::vector<int> v,const int i) // v is assumed to be sorted
{
	// Increase size by one
	v.resize(v.size() + 1);

	// Case i > last of v
	if (v[v.size()-2] < i){
		v[v.size()-1] = i;
		return v;
	}
	
	v[v.size()-1] = v[v.size()-2] + 1; // So v is still sorted

	// Find where to place the new value
	auto upper = std::upper_bound(v.begin(), v.end(), i);
	
	// Shift the elements above i to the right
	for(std::vector<int>::reverse_iterator rit = v.rbegin(); &*rit > &*upper ; ++rit)
		*rit = *(rit+1);

	// Insert new value
	
	*upper = i;

	return v;
}

//(Vector?)
void BronKerbosch::unite_vector_element_ref(std::vector<int>& v,const int i) // v is assumed to be sorted
{
	// Increase size by one
	v.resize(v.size() + 1);

	// Case i > last of v
	if (v[v.size()-2] < i){
		v[v.size()-1] = i;
		return;
	}

	v[v.size()-1] = v[v.size()-2] + 1; // So v is still sorted

	// Find where to place the new value
	auto upper = std::upper_bound(v.begin(), v.end(), i);

	// Shift the elements above i to the right
	for(std::vector<int>::reverse_iterator rit = v.rbegin(); &*rit > &*upper ; ++rit)
		*rit = *(rit+1);

	// Insert new value
	*upper = i;	
}


void BronKerbosch::find_cliques_withoutPivot(std::vector<int> clique, std::vector<int> candidates, std::vector<int> banned)
{
	bronkerbosch_calls_display(clique, candidates, banned);	

	if(candidates.empty() && banned.empty())
		add_max_clique(clique);
	
	// TREAT CANDIDATES BACKWARD
	for(std::vector<int>::const_reverse_iterator rit = candidates.crbegin(); rit != candidates.crend(); ++rit) {
		// RECURSIVE CALL
		find_cliques_withoutPivot(unite_vector_element_val(clique, *rit), intersect(candidates, g.adjacency_OLD[*rit]), intersect(banned, g.adjacency_OLD[*rit]));

		// P := \ {v}
		candidates.pop_back(); // MUCH BETTER TO TREAT CANDIDATES BACKWARDS AND USE POP_BACK !

		// X := X U {v} // Need to be sorted -> Consuming 
		banned.push_back(*rit); 
		std::sort (banned.begin(), banned.end());

		std::cout << "candidates left: ";
		for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
			std::cout << " " << *it2+1;}
		std::cout << std::endl;

		std::cout << "banned : ";
		for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3+1;}
		std::cout << std::endl;
	}		
}

// Find the vertice with max degree within two vectors of vertices
// Shouldn't we record degrees at the begining and access them without search ? 
int BronKerbosch::find_pivot(std::vector<int> const& v1, std::vector<int> const& v2)
{	
	int degree;
	int max_degree = -1;
	int pivot = -1;
	for(std::vector<int>::const_iterator it = v1.begin(); it != v1.end(); ++it){
		degree = g.adjacency_OLD[*it].size();		
		if (degree > max_degree){ max_degree = degree; pivot = *it;}}
	for(std::vector<int>::const_iterator it = v2.begin(); it != v2.end(); ++it){
		degree = g.adjacency_OLD[*it].size();
		if (degree > max_degree){ max_degree = degree; pivot = *it;}}
	return pivot;		
}

void BronKerbosch::find_cliques_withPivot(std::vector<int> clique, std::vector<int> candidates, std::vector<int> banned)
{
	bronkerbosch_calls_display(clique, candidates, banned);	

	if(candidates.empty() && banned.empty())
		add_max_clique(clique);

	// Choose a pivot among P U X with max_degree
	int pivot = find_pivot(candidates, banned);
	std::cout << "Pivot : " << pivot+1 << std::endl;
	
	std::cout << " Candidates - {";
	for(std::vector<int>::const_iterator it = g.adjacency_OLD[pivot].begin(); it != g.adjacency_OLD[pivot].end(); ++it){
		std::cout << *it+1 << " ";
	}
	std::cout <<"}"<< std::endl;

	// Restrict candidates to P \ N(pivot)
	std::vector<int> restricted_candidates;	
	restricted_candidates.resize(candidates.size()); // Necessary ?	
	restricted_candidates = eliminate(candidates, g.adjacency_OLD[pivot]);
	std::cout << "Ok " << std::endl;

	std::cout << "Restricted candidates in P - N("<<pivot+1<<") : {";
	for(std::vector<int>::const_iterator it = restricted_candidates.begin(); it != restricted_candidates.end(); ++it){
		std::cout << " " << *it+1;}
	std::cout <<"}"<< std::endl;

	// TREAT CANDIDATES BACKWARD
	for(std::vector<int>::const_reverse_iterator rit = restricted_candidates.crbegin(); rit != restricted_candidates.crend(); ++rit) {
		// RECURSIVE CALL
		find_cliques_withPivot(unite_vector_element_val(clique, *rit), intersect(candidates, g.adjacency_OLD[*rit]), intersect(banned, g.adjacency_OLD[*rit]));
		// P := \ {v}
		candidates.pop_back(); // MUCH BETTER TO TREAT CANDIDATES BACKWARDS AND USE POP_BACK !

		// X := X U {v} // Need to be sorted -> Consuming 
		banned.push_back(*rit); 
		std::sort (banned.begin(), banned.end());

		std::cout << "candidates left: ";
		for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
			std::cout << " " << *it2+1;}
		std::cout << std::endl;

		std::cout << "banned : ";
		for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3+1;}
		std::cout << std::endl;
	}		
}

void BronKerbosch::order_by_degree()
{	
	by_degree.resize(g.capacity());
	degree.resize(g.capacity());
	for (auto v : g.nodes) {
        degree[v] = g.adjacency_OLD[v].size();
        //by_degree[degree[v]].push_back(v); FALSE
    }

	for(std::vector<int>::const_iterator it = by_degree.begin(); it != by_degree.end(); ++it){
	std::cout << " " << *it;}
	std::cout << std::endl;
}


/**************************************/
/*			Display methods 		  */
/**************************************/
void BronKerbosch::display()
{
	std::cout << " Maximal cliques " << std::endl;
	int i = 1;
	for(std::vector<std::vector<int>>::const_iterator it = maximal_cliques.begin(); it != maximal_cliques.end(); ++it) {
		if (!(*it).empty()){
			std::cout << " Clique " << i << " { ";
			i++;
			for (std::vector<int>::const_iterator itN = (*it).begin(); itN < (*it).end(); ++itN) {
				std::cout << *itN+1 << " ";
			}
			std::cout << " } size " << (*it).size() <<std::endl;
		}
	}
}

void BronKerbosch::bronkerbosch_calls_display(std::vector<int> const& clique, std::vector<int> const& candidates, std::vector<int> const& banned)
{
	// Display calls
	std::cout << "BronKerbosch({" ;
	for(std::vector<int>::const_iterator it1 = clique.begin(); it1 != clique.end(); ++it1){
		std::cout << " " << *it1+1;}
	std::cout << " }, {";
	for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
		std::cout << " " << *it2+1;}
	std::cout << " }, {";
	for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3+1;}
	std::cout << " })"<< std::endl;
}

// Add the print_container



// degeneracy_finder::degeneracy_finder(const graph& g)
//     : g(g)
// {
//
//      degrees.resize(g.capacity());
//      iterators.resize(g.capacity());
//      ordered.initialise(0, g.capacity()-1, bitset::empt);
// }
//
//
//
// void degeneracy_finder::get_degeneracy_order( std::vector< int >& order  ) {
//              for (auto v : g.nodes) {
//                  // auto vd = g.neighbor(v).size();
//                  if (vd >= buckets.size())
//                      buckets.resize(vd + 1);
//                  buckets[vd].push_front(v);
//                  degrees[v] = vd;
//                  iterators[v] = buckets[vd].begin();
//              }
//
//              ordered.clear();
//     while (true) {
//         size_t i{0};
//         for (; i != buckets.size(); ++i)
//             if (!buckets[i].empty())
//                 break;
//         if (i == buckets.size())
//             break;
//         auto v = buckets[i].back();
//         order.push_back(v);
//         buckets[i].pop_back();
//         ordered.fast_add(v);
//         for (auto u : neighbor(v)) {
//             if (ordered.fast_contain(u))
//                 continue;
//             auto& ud = degrees[u];
//             buckets[ud].erase(iterators[u]);
//             --ud;
//             buckets[ud].push_front(u);
//             iterators[u] = buckets[ud].begin();
//         }
//     }
// }


} // namespace gc

