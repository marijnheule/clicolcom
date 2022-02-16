#ifndef __CG_GRAPH_HH
#define __CG_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"
#include "utils.hpp"

#include <minicsp/core/utils.hpp>
#include <minicsp/mtl/Heap.h>

#include <algorithm>
#include <list>
#include <random>
#include <vector>

// #define _DEBUG_CLIQUE
// #define SUBCLASS

// #define _DEBUG_PRUNING
// #define _FULL_PRUNING

namespace gc
{

using weight = int64_t;
using bitset = BitSet;
using edge = std::pair<int, int>;

struct arc {
    int v[2];
    const int operator[](const int i) { return v[i]; }
};

template <class adjacency_struct> class basic_graph
{
public:
    bitset nodeset;
    intstack nodes;
    std::vector<adjacency_struct> matrix;

    basic_graph() {}
    explicit basic_graph(int nv)
        : nodeset(0, nv - 1, bitset::full)
        // : nodeset(0, nv - 1, bitset::empt)
        , matrix(nv)
    {
        nodes.reserve(nv);
        nodes.fill();
        for (auto& bs : matrix) {
            bs.initialise(0, nv, bitset::empt);
        }
    }

    basic_graph& operator=(const basic_graph&) = default;

    template <typename other_struct>
    typename std::enable_if<!std::is_same<other_struct, adjacency_struct>{},
        basic_graph>::type&
    operator=(const basic_graph<other_struct>& g)
    {
        nodes.reserve(g.capacity());
        nodeset.initialise(0, g.capacity() - 1, bitset::empt);
        matrix.resize(g.capacity());
        for (auto v : g.nodes) {
            matrix[v].initialise(0, g.capacity() - 1, bitset::empt);
            add_node(v);
            for (auto u : g.matrix[v]) {
                matrix[v].add(u);
            }
        }
        return *this;
    }

    basic_graph(const basic_graph&) = default;
    template <typename other_struct,
        typename E = typename std::enable_if<
            !std::is_same<other_struct, adjacency_struct>{}>::type>
    basic_graph(const basic_graph<other_struct&> g)
    {
        this->operator=(g);
    }
    template <class other_struct>
    basic_graph(const basic_graph<other_struct>& g, std::vector<int>& vmap)
    {
        if (g.size() > 0) {
            nodes.reserve(g.size());
            nodes.fill();
            nodeset.initialise(0, g.size() - 1, bitset::full);
            matrix.resize(g.size());
            int i = 0;
            for (auto v : g.nodes) {
                matrix[i].initialise(0, g.size() - 1, bitset::empt);
                vmap[v] = i++;
            }
            assert(i == g.size());
            for (auto v : g.nodes) {
                for (auto u : g.matrix[v]) {
                    if (g.nodeset.fast_contain(u))
                        matrix[vmap[v]].add(vmap[u]);
                }
            }
            canonize();
        }
    }
    basic_graph(basic_graph&&) = default;
    basic_graph& operator=(basic_graph&&) = default;

    // Method to initialize a non-const basic graph from a const basic_graph

    int capacity() const { return matrix.size(); }
    int size() const { return nodes.size(); }
    long count_edges() const;
    bool has_edge(const int u, const int v) const;

    void add_edge(const int u, const int v);

    void add_edges(const int u, const adjacency_struct& N);

    void add_node(const int v);

    void add_clique(const adjacency_struct& C);

    void remove_edge(int u, int v);

    void remove_node(int v);

    void remove(const adjacency_struct& nodes);

    void remove(const int v);

    void clear();

    void canonize();

    void describe(std::ostream& os, const int num_edges = -1) const;

    // whether v is simplicial (if not, a witness edge is cpied)
    bool simplicial(
        const int v, std::pair<int, int>& witness, bitset& util_set) const;

    void check_basic_consistency() const
    {
        std::cout << "checking basic-graph consistency" << std::endl;

        assert(nodes.size() == nodeset.size());
        for (auto v : nodes) {
            assert(!(matrix[v].fast_contain(v)));

            assert(nodeset.fast_contain(v));
        }

        for (auto v : nodes) {
            for (auto u : matrix[v]) {

                if (!matrix[u].fast_contain(v)) {
                    std::cout << v << " has a non-symmetric neighbor (" << u
                              << ")\n"
                              << "N(" << v << ") = " << matrix[v] << "\nN(" << u
                              << ") = " << matrix[u] << "\n";
                }

                if (!nodeset.fast_contain(u)) {
                    std::cout << v
                              << " has a neighbor that is not in the nodeset ("
                              << u << ")\n"
                              << nodeset << std::endl;
                }

                assert(matrix[u].fast_contain(v));
                assert(nodeset.fast_contain(u));
            }
        }
    }
};

template <class adjacency_struct>
class graph : public basic_graph<adjacency_struct>
{
public:
    using basic_graph<adjacency_struct>::matrix;
    using basic_graph<adjacency_struct>::nodes;
    using basic_graph<adjacency_struct>::nodeset;

    // we keep a copy of the original matrix because we modify matrix
    // when we do merge/separate
    std::vector<adjacency_struct> origmatrix;

    // checkpointing
    int cur_ckpt{0};
    std::vector<int> removed;

    // trail of removed edges, including a lim which delimits
    // checkpoints
    std::vector<std::pair<int, int>> edge_trail;
    std::vector<std::size_t> edge_lim;

    // the partitions generated by merging vertices. initially each
    // partition is a singleton, so rep_of[i]=i, partition[i] =
    // {i}. When we merge u and v, we keep one vertex (say u) as
    // representative of the partition, set partition[u] to be the
    // union of partition[u] and partition[v] (which were previously
    // disjoint) and set rep_of[i]=u for all i in partition[v],
    // including v. At the end g has u in its vertex set, but none of
    // the other vertices of partition[u].
    std::vector<int> rep_of;
    std::vector<std::vector<int>> partition;
    // information to backtrack the above: previous rep_of and
    // previous partition sizes
    std::vector<std::vector<size_t>> partition_size_trail;
    std::vector<std::vector<int>> rep_of_trail;

    // buffers
    bitset util_set, diff2;
    bitset partu, partv;

    //--------------------------------------------------
    // private, but out in the open

    // add an edge that will be backtracked later
    void add_dirty_edge(int u, int v);

public:
    graph()
        : basic_graph<adjacency_struct>()
    {
    }
    explicit graph(int nv)
        : basic_graph<adjacency_struct>(nv)
        , origmatrix(nv)
        , rep_of(nv)
        , partition(nv)
        , util_set(0, nv - 1, bitset::empt)
        , diff2(0, nv - 1, bitset::empt)
        , partu(0, nv - 1, bitset::empt)
        , partv(0, nv - 1, bitset::empt)
    {
        for (auto& bs : origmatrix) {
            bs.initialise(0, nv, bitset::empt);
        }
        for (auto v : nodes) {
            rep_of[v] = v;
            partition[v].push_back(v);
        }
    }

    void init_structures()
    {
        origmatrix.resize(this->capacity());
        rep_of.resize(this->capacity());
        partition.resize(this->capacity());
        util_set.initialise(0, this->capacity() - 1, bitset::empt);
        diff2.initialise(0, this->capacity() - 1, bitset::empt);
        partu.initialise(0, this->capacity() - 1, bitset::empt);
        partv.initialise(0, this->capacity() - 1, bitset::empt);

        for (auto v : nodes) {
            origmatrix[v].initialise(0, this->capacity() - 1, bitset::empt);
            origmatrix[v].copy(matrix[v]);
            rep_of[v] = v;
            partition[v].push_back(v);
        }
    }
    graph& operator=(const graph&) = default;
    template <typename other_struct>
    typename std::enable_if<!std::is_same<other_struct, adjacency_struct>{},
        graph>::type&
    operator=(const graph<other_struct>& g)
    {
        this->basic_graph<adjacency_struct>::operator=(g);
        init_structures();

        return *this;
    }

    graph(const graph&) = default;
    template <typename other_struct,
        typename E = typename std::enable_if<
            !std::is_same<other_struct, adjacency_struct>{}>::type>
    graph(const graph<other_struct>& g)
    //: basic_graph<adjacency_struct>(g)
    {
        this->operator=(g);
    }
    template <class other_struct>
    graph(const graph<other_struct>& g, std::vector<int>& vmap)
        : basic_graph<adjacency_struct>(g, vmap)
    {
        init_structures();
    }
    graph(graph&&) = default;
    graph& operator=(graph&&) = default;

    int capacity() const { return matrix.size(); }

    void add_edge(int u, int v);

    // merge vertices and return the id of the new vertex (one of
    // u, v)
    int merge(int u, int v);

    // separate u and v. Just adds an edge, but it is reversible
    // through checkpointing
    void separate(int u, int v);

    int contractPreprocess();

    int checkpoint();

    void restore(int ckpt);

    int current_checkpoint() const { return cur_ckpt; }

    int representative_of(const int v) const { return rep_of[v]; }

    // debugging
    void tell_class() const { std::cout << "GRAPH\n"; }

    void check_consistency() const;
};

template <class adjacency_struct> struct clique_finder {
    const graph<adjacency_struct>& g;
    std::vector<adjacency_struct> cliques;
    std::vector<int> clique_sz;
    // std::vector<adjacency_struct> candidates;
    std::vector<gc::bitset> candidates;
    std::vector<int> last_clique;
    int num_cliques;
    int limit;

    std::vector<int> util_vec;

    // for pruning
    bool update_nn;

#ifdef _DEBUG_PRUNING
    std::vector<std::vector<int>> num_neighbors_of;
    std::vector<int> pruning;
#endif

    std::vector<int> new_tight_cliques;
    std::vector<std::vector<int>> neighbor_count;
    std::vector<int> col_pruning;

#ifdef _FULL_PRUNING
    std::vector<std::vector<int>> cliques_of;
#endif

    clique_finder(const graph<adjacency_struct>& ig, const int c = 0xfffffff)
        : g(ig)
        , num_cliques(1)
        , limit(c)
        , update_nn(false)
    {
        auto m = std::min(limit, g.capacity());
        last_clique.resize(g.capacity());
        cliques.resize(m);
        clique_sz.resize(m);
        candidates.resize(m);
        for (auto& b : cliques)
            b.initialise(0, g.capacity(), bitset::empt);
        for (auto& b : candidates)
            b.initialise(0, g.capacity(), bitset::full);
    }

    void switch_pruning()
    {
#ifdef _FULL_PRUNING
        cliques_of.resize(g.capacity());
#endif
        update_nn = true;
    }

#ifdef _DEBUG_PRUNING
    void compute_neighbors_in(const int i, const int sz)
    {
        for (auto v : cliques[i]) {
            for (auto u : g.matrix[v]) {
                if (g.nodeset.fast_contain(u) and !cliques[i].fast_contain(u)) {

                    if (++num_neighbors_of[i][u] == sz) {
                        pruning.push_back(i);
                        pruning.push_back(u);
                    }
                }
            }
        }
    }
#endif

    void count_neighbors(const int i)
    {
        for (auto v : cliques[i]) {
            for (auto u : g.matrix[v]) {
                if (g.nodeset.fast_contain(u) and !cliques[i].fast_contain(u)) {
                    if (++neighbor_count[i][u] == clique_sz[i] - 1) {
                        col_pruning.push_back(i);
                        col_pruning.push_back(u);
                    }
                }
            }
        }
    }

#ifdef _DEBUG_PRUNING
    void notify_ub(const int ub)
    {
        pruning.clear();

        int n
            = std::min(static_cast<int>(num_neighbors_of.size()), num_cliques);
        for (int i = 0; i < n; ++i) {
            std::fill(begin(num_neighbors_of[i]), end(num_neighbors_of[i]), 0);
        }
        if (n < num_cliques) {
            num_neighbors_of.resize(num_cliques);
            for (int i = n; i < num_cliques; ++i) {
                num_neighbors_of[i].resize(g.capacity(), 0);
            }
        }

        for (auto i{0}; i < num_cliques; ++i) {
            if (clique_sz[i] == ub - 1) {
                compute_neighbors_in(i, ub - 2);
            }
        }
    }
#endif

    void compute_pruning(const int maxcs, std::vector<arc>& changed_edges,
        std::vector<int>& changed_vertices)
    {
        col_pruning.clear();

        // compute the count for new cliques
        for (auto cl : new_tight_cliques) {
            assert(cl < num_cliques);
            if (cl >= static_cast<int>(neighbor_count.size())) {
                neighbor_count.resize(cl + 1);
            }
            if (neighbor_count[cl].size() == 0)
                neighbor_count[cl].resize(g.capacity(), 0);
            else
                std::fill(
                    begin(neighbor_count[cl]), end(neighbor_count[cl]), 0);

            count_neighbors(cl);
        }

// new edges with old cliques
#ifdef _FULL_PRUNING
        if (new_tight_cliques.size() < num_cliques) {
            for (auto v : changed_vertices) {
                cliques_of[v].clear();
            }

            int j = 0;
            for (int i = 0; i < num_cliques; ++i) {
                if (new_tight_cliques.size() > j
                    and i == new_tight_cliques[j]) {
                    ++j;
                } else {
                    if (clique_sz[i] == maxcs) {
                        for (auto v : changed_vertices)
                            if (cliques[i].fast_contain(v))
                                cliques_of[v].push_back(i);
                    }
                }
            }

            for (auto e : changed_edges) {
                for (int i = 0; i < 2; ++i) {
                    for (auto cl : cliques_of[e[1 - i]]) {
                        if (!cliques[cl].fast_contain(e[i])
                            and ++neighbor_count[cl][e[i]] == maxcs - 1) {
                            col_pruning.push_back(cl);
                            col_pruning.push_back(e[i]);
                        }
                    }
                }
            }
        }
#endif

        // new_tight_cliques.clear();
    }

    // clear previously cached results
    void clear()
    {
        num_cliques = 0;
    }
    // initialize a new clique
    void new_clique();

    // initialize a new color
    void new_color();

    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);

    // insert v into the col^th color. assumes it fits. Puts vertices
    // added from candidates[i] into diff
    void insert_color(int v, int clq, bitset& diff);

    // heuristically find a set of cliques and return the size of the
    // largest

    template <class ordering>
    int find_cliques(
        ordering o, const int lower, const int upper, const int l = 0xfffffff)
    {
        // new_tight_cliques.clear();
        if (l < limit)
            limit = l;

        clear();
        if (o.size() == 0)
            return 0;
        for (auto u : o) {
            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(u)) {
                    found = true;
                    insert(u, i);
                }
            if (!found && num_cliques < limit) {
                new_clique();
                insert(u, num_cliques - 1);
            }
        }

        for (auto u : o) {
            for (int i = last_clique[u] + 1; i < num_cliques; ++i)
                if (candidates[i].fast_contain(u)) {
                    insert(u, i);
                }
        }

        auto maxclique{*std::max_element(
            begin(clique_sz), begin(clique_sz) + num_cliques)};

        new_tight_cliques.clear();

        if (update_nn and lower < upper - 1 and maxclique == upper - 1) {
            for (auto i{0}; i < num_cliques; ++i) {
                if (clique_sz[i] == maxclique) {
                    new_tight_cliques.push_back(i);

                    // std::cout << " " << i;
                }
            }
            // if(new_tight_cliques.size() > 0)
            // 		std::cout << std::endl;
        }

        return maxclique;
    }

    void filter_cliques(int cutoff, int tight)
    {

        // when we filter, we must update all fields
        int i{0}, j{0}, k{0};
        for (; i != num_cliques; ++i)
            if (clique_sz[i] >= cutoff) {

                // rename the new cliques in case of filtering
                if (update_nn and clique_sz[i] == tight
                    and static_cast<int>(new_tight_cliques.size()) > k
                    and new_tight_cliques[k] == i) {
                    new_tight_cliques[k] = j;
                    ++k;
                }

                using std::swap;
                swap(cliques[i], cliques[j]);
                swap(clique_sz[i], clique_sz[j]);
                swap(candidates[i], candidates[j]);

                ++j;
            }
        num_cliques = j;
    }

    void remap_clique_to_reps(bitset& clq, bitset& N)
    {
        util_vec.clear();
        std::copy(begin(clq), end(clq), back_inserter(util_vec));
        N.clear();
        bool first{true};
        for (auto u : util_vec) {
            if (g.rep_of[u] != u) {
                clq.fast_remove(u);
                clq.fast_add(g.rep_of[u]);
            }
            if (first) {
                N.copy(g.matrix[g.rep_of[u]]);
                first = false;
            } else {
                N.intersect_with(g.matrix[g.rep_of[u]]);
            }
        }
    }

    int extend_cliques(
        const std::vector<int>& cand, const int lower, const int upper)
    {
        int lb = -1;
        for (int i = 0; i != num_cliques; ++i) {
            remap_clique_to_reps(cliques[i], candidates[i]);

            auto sz_before{clique_sz[i]};

            clique_sz[i] = extend_clique(cliques[i], candidates[i], cand);

            if (update_nn and sz_before < clique_sz[i] and clique_sz[i] == upper - 1) {
                new_tight_cliques.push_back(i);
            }

            lb = std::max(lb, clique_sz[i]);
            if (lb >= upper)
                break;
        }

        return lb;
    }

    int extend_clique(bitset& clq, bitset& N, const std::vector<int>& cand)
    {
        for (auto u : cand) {
            if (N.fast_contain(u)) {
                clq.fast_add(u);
                N.intersect_with(g.matrix[u]);
            }
        }
        return clq.size();
    }

    // template <class ordering> int find_cliques(std::vector<int>::iterator& b,
    // std::vector<int>::iterator& e, const int l=0xfffffff)
    // {
    //     if (l < limit)
    //         limit = l;
    //     clear();
    //     if (e == b)
    //         return 0;
    //
    //     for (auto ui=b; ui!=e; ++ui) {
    //         bool found{false};
    //                                          int u = *ui;
    //         for (int i = 0; i != num_cliques; ++i)
    //             if (candidates[i].fast_contain(u)) {
    //                 found = true;
    //                 insert(u, i);
    //             }
    //         if (!found && num_cliques < limit) {
    //             new_clique();
    //             insert(u, num_cliques - 1);
    //         }
    //     }
    //
    //     for (auto ui=b; ui!=e; ++ui) {
    //                                          int u = *ui;
    //         for (int i = last_clique[u] + 1; i < num_cliques; ++i)
    //             if (candidates[i].fast_contain(u)) {
    //                 insert(u, i);
    //             }
    //     }
    //
    //     return *std::max_element(
    //         begin(clique_sz), begin(clique_sz) + num_cliques);
    // }

    template <class ordering>
    void find_clique_cover(ordering o)
    {
        clear();
        if (o.size() == 0)
            return;
        for (auto u : o) {
            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(u)) {
                    found = true;
                    insert(u, i);
                    break;
                }
            if (!found) {
                new_clique();
                insert(u, num_cliques - 1);
            }
        }
    }

    void sort_cliques(const int size)
    {
        for (auto cl = 0; cl < num_cliques; ++cl)
            if (clique_sz[cl] >= size)
                cliques[cl].canonize();
    }
};

template <class graph_struct> struct degeneracy_finder {

    const graph_struct& g;
    int degeneracy;
    std::vector<int> order;
    std::vector<int> degrees;
    std::vector<std::list<int>::iterator> iterators;
    std::vector<bool> ordered;
    std::vector<std::list<int>> buckets;

    std::vector<std::vector<int>::iterator> core;
    std::vector<int> core_degree;

    degeneracy_finder(const graph_struct& g)
        : g(g)
        , degeneracy{-1}
        , degrees(g.capacity())
        , iterators(g.capacity())
        , ordered(g.capacity())
    {
    }

    void degeneracy_ordering(); // Matula & Beck
    void clear();
    void display_ordering();
};

template <class adjacency_struct> struct minimal_triangulator {
    const basic_graph<adjacency_struct>& g;
    basic_graph<adjacency_struct> g_filled;
    basic_graph<adjacency_struct> g_section;
    std::vector<std::pair<int, int>> fill_edges;

    std::vector<int> increasing_ordering;
    std::vector<int> minimum_degree_ordering;

    minimal_triangulator(const basic_graph<adjacency_struct>& g)
        : g(g)
        , g_filled(g.capacity())
        , g_section(g.capacity())
    {
        // COPY G
        //      g_filled = g;
        //      g_section = g;
        for (auto v : g.nodes) {
            g_filled.add_node(v);
            g_filled.matrix[v].copy(g.matrix[v]);

            g_section.add_node(v);
            g_section.matrix[v].copy(g.matrix[v]);
        }

        for (auto u : g.nodes) {
            increasing_ordering.push_back(u);
        }
    }
    // Sorting nodes by degree
    void stable_sort();

    void section_node(int n);

    void elimination_game(const std::vector<int>& ordering);
};

/** IMPLEMENTATION **/

// template<class graph_struct>

template <class adjacency_struct>
long basic_graph<adjacency_struct>::count_edges() const
{
    long m{0};
    for (auto v : nodes)
        m += matrix[v].size();
    return m / 2;
}

template <class adjacency_struct>
bool basic_graph<adjacency_struct>::has_edge(const int u, const int v) const
{
    bool he{matrix[u].fast_contain(v)};
    if (he)
        assert(matrix[v].fast_contain(u));
    return he;
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::add_edge(const int u, const int v)
{
    matrix[u].add(v);
    matrix[v].add(u);
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::add_edges(
    const int u, const adjacency_struct& N)
{
    matrix[u].union_with(N);
    for (auto v : N) {
        matrix[v].add(u);
    }
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::add_node(const int v)
{
    nodes.add(v);
    nodeset.add(v);
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::add_clique(const adjacency_struct& C)
{
    for (auto v : C) {
        assert(!nodeset.fast_contain(v));

        add_node(v);
        matrix[v].union_with(C);
        matrix[v].remove(v);
    }
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::remove(const adjacency_struct& toremove) {
    for (auto v : toremove) {
        nodes.remove(v);
    }
    nodeset.setminus_with(toremove);
    if (nodes.size() < toremove.size()) {
        for (auto v : nodes) {
            matrix[v].intersect_with(nodeset);
        }
    } else {
        for (auto v : nodes)
            matrix[v].setminus_with(toremove);
    }
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::remove(const int v)
{
    remove_node(v);

    for (auto u : matrix[v])
        matrix[u].remove(v);
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::remove_edge(int u, int v)
{
    matrix[u].remove(v);
    matrix[v].remove(u);
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::remove_node(int v)
{
    nodes.remove(v);
    nodeset.fast_remove(v);
}

template <class adjacency_struct> void basic_graph<adjacency_struct>::clear()
{
    for (auto v : nodes) {
        matrix[v].clear();
    }
    nodeset.clear();
    nodes.clear();
}

template <class adjacency_struct> void basic_graph<adjacency_struct>::canonize()
{
    for (auto v : nodes) {
        matrix[v].canonize();
        //         matrix[v].sort();
        // std::unique(begin(matrix[v]), end(matrix[v]));
    }
}

template <class adjacency_struct>
bool basic_graph<adjacency_struct>::simplicial(
    const int v, std::pair<int, int>& witness, bitset& util_set) const
{
    for (auto u : matrix[v]) {
        util_set.copy(matrix[v]);
        util_set.setminus_with(matrix[u]);
        if (util_set.empty())
            return true;
        witness.first = u;
        witness.second = util_set.min();
    }
    return false;
}

template <class adjacency_struct>
void graph<adjacency_struct>::add_dirty_edge(int u, int v)
{
    if (matrix[u].fast_contain(v))
        return;
    for (auto vp : partition[v]) {
        matrix[u].fast_add(vp);
        if (cur_ckpt > 0)
            edge_trail.push_back({u, vp});
    }
}

template <class adjacency_struct>
void graph<adjacency_struct>::add_edge(int u, int v)
{

    if (matrix.size() <= static_cast<size_t>(v)) {
        std::cout << "ERROR: try to add edge " << u << "-" << v << " / "
                  << nodes.size() << " " << nodeset.size() << " "
                  << matrix.size() << std::endl;
    }

    assert(matrix.size() > static_cast<size_t>(v));
    assert(matrix.size() > static_cast<size_t>(u));
    assert(origmatrix.size() > static_cast<size_t>(v));
    assert(origmatrix.size() > static_cast<size_t>(u));

    matrix[u].add(v);
    matrix[v].add(u);
    origmatrix[u].add(v);
    origmatrix[v].add(u);
}

template <class adjacency_struct>
int graph<adjacency_struct>::merge(int u, int v)
{
    // util_set.clear();
    util_set.copy(matrix[v]);
    util_set.setminus_with(matrix[u]);

    // diff2.clear();
    diff2.copy(matrix[u]);
    diff2.setminus_with(matrix[v]);

    matrix[u].union_with(matrix[v]);
    for (auto w : util_set) {
        add_dirty_edge(w, u);
        add_dirty_edge(u, w);
    }
    for (auto w : diff2) {
        for (auto vp : partition[v])
            add_dirty_edge(w, vp);
    }

    // update rep_of for the partition that was absorbed
    for (auto vp : partition[v])
        rep_of[vp] = u;
    std::copy(
        begin(partition[v]), end(partition[v]), back_inserter(partition[u]));

    nodes.remove(v);
    nodeset.fast_remove(v);
    removed.push_back(v);

    return u;
}

// separate u and v. Just adds an edge, but it is reversible through
// checkpointing
template <class adjacency_struct>
void graph<adjacency_struct>::separate(int u, int v)
{
    if (matrix[u].contain(v)) {
        assert(matrix[v].contain(u));
        return;
    }
    // std::cout << "separating " << u << " from " << v << std::endl;
    add_dirty_edge(u, v);
    add_dirty_edge(v, u);
}

template <class adjacency_struct>
int graph<adjacency_struct>::contractPreprocess()
{
    int num_contractions = 0;
    bool some_propagation = true;
    while (some_propagation) {
        some_propagation = false;
        for (auto u : nodes) {
            for (auto v : nodes) {
                if (u != v && !matrix[u].fast_contain(v)) {
                    util_set.copy(matrix[v]);
                    util_set.setminus_with(matrix[u]);
                    if (!util_set.intersect(nodeset)) {
                        // N(v) <= N(U)s
                        some_propagation = true;
                        ++num_contractions;
                        merge(u, v);
                    }
                }
            }
        }
    }
    return num_contractions;
}

template <class adjacency_struct> int graph<adjacency_struct>::checkpoint()
{
    ++cur_ckpt;

    if (static_cast<size_t>(cur_ckpt) >= rep_of_trail.size()) {
        // make space for the copying part
        rep_of_trail.resize(cur_ckpt);
        partition_size_trail.resize(cur_ckpt);
        partition_size_trail[cur_ckpt - 1].resize(capacity());
    }
    // trailing
    edge_lim.push_back(edge_trail.size());

    // copying
    rep_of_trail[cur_ckpt - 1] = rep_of;
    for (auto v : nodes)
        partition_size_trail[cur_ckpt - 1][v] = partition[v].size();

    return cur_ckpt;
}

template <class adjacency_struct>
void graph<adjacency_struct>::restore(int ckpt)
{
    cur_ckpt = ckpt;

    int trailnewsize = edge_lim[ckpt];
    for (int i = trailnewsize; static_cast<unsigned>(i) != edge_trail.size();
         ++i) {
        auto e = edge_trail[i];
        matrix[e.first].fast_remove(e.second);
        matrix[e.second].fast_remove(e.first);
    }
    edge_trail.resize(trailnewsize);
    edge_lim.resize(ckpt);

    rep_of = rep_of_trail[cur_ckpt];

    while (!removed.empty()) {
        auto v = removed.back();
        assert(!nodes.contain(v));
        if (rep_of[v] != v)
            break;
        nodes.add(v);
        nodeset.fast_add(v);
        removed.pop_back();
    }

    for (auto v : nodes)
        partition[v].resize(partition_size_trail[cur_ckpt][v]);
}

template <class adjacency_struct>
void basic_graph<adjacency_struct>::describe(std::ostream& os, const int num_edges) const
{
    int m = num_edges;
    if (m < 0) {
        m = count_edges();
    }

    os << "#vertices = " << this->size() << ",  #edges = " << m
       << ",  density = "
       << (m > 0
                  ? (double)(2 * m)
                      / ((double)(this->size()) * (double)(this->size() - 1))
                  : 1)
            * 100
       << "%";
}

template <class adjacency_struct>
void graph<adjacency_struct>::check_consistency() const
{
    std::cout << "checking graph consistency" << std::endl;

    // basic_graph<adjacency_struct>::check_consistency();

    for (int i = 0; i != capacity(); ++i)
        assert((!nodes.contain(i) or nodeset.contain(i))
            and (!nodeset.contain(i) or nodes.contain(i)));

    for (int i = 0; i != capacity(); ++i)
        assert((!nodes.contain(i) || rep_of[i] == i)
            && (rep_of[i] != i || nodes.contain(i)));

    bitset bs(0, capacity() - 1, bitset::full);
    for (auto v : removed) {
        assert(!nodes.contain(v));
        bs.fast_remove(v);
    }
    assert(bs == nodeset);

    for (auto v : nodes) {
        assert(!partition[v].empty());
        assert(partition[v][0] == v);
        assert(rep_of[v] == v);
        bs.clear();
        for (auto u : partition[v]) {
            assert(rep_of[u] == v);
            bs.union_with(origmatrix[u]);
        }
        if (!bs.included(matrix[v])) {
            std::cout << "    bs[" << v << "] = " << bs << "\n";
            std::cout << "matrix[" << v << "] = " << matrix[v] << "\n";
        }
        assert(bs.included(matrix[v]));
    }

    for (auto v : nodes) {
        for (auto u : nodes) {
            if (u == v)
                continue;
            // if (u < v)
            //     continue;
            if (matrix[v].fast_contain(u)) {
                assert(matrix[u].fast_contain(v));
                for (auto vp : partition[v])
                    if (!matrix[u].fast_contain(vp)) {
                        std::cout << "u = " << u << " v = " << v
                                  << " vp = " << vp << "\n";
                        std::cout << "matrix[v] = " << matrix[v] << std::endl;
                        std::cout << "matrix[u] = " << matrix[u] << std::endl;
                        std::cout
                            << "partition[v] = "
                            << print_container<std::vector<int>>{partition[v]}
                            << std::endl;
                        std::cout
                            << "partition[u] = "
                            << print_container<std::vector<int>>{partition[u]}
                            << std::endl;
                        assert(matrix[u].fast_contain(vp));
                    }
                for (auto up : partition[u]) {
                    if (!matrix[v].fast_contain(up)) {
                        std::cout << "matrix[v] = " << matrix[v] << std::endl;
                        std::cout << "matrix[u] = " << matrix[u] << std::endl;
                        assert(matrix[v].fast_contain(up));
                    }
                }
            }
        }
    }
}

template <class adjacency_struct>
void clique_finder<adjacency_struct>::new_clique()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].copy(g.nodeset);
    ++num_cliques;

    // if (update_nn and num_neighbors_of.size() < num_cliques) {
    //     num_neighbors_of.resize(num_cliques);
    //     num_neighbors_of.back().resize(g.capacity(), 0);
    // }
}
// initialize a new color
template <class adjacency_struct>
void clique_finder<adjacency_struct>::new_color()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].clear();
    ++num_cliques;
}
// insert v into the clq^th clique. assumes it fits
template <class adjacency_struct>
void clique_finder<adjacency_struct>::insert(int v, int clq)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    candidates[clq].intersect_with(g.matrix[v]);
    last_clique[v] = clq;

    // if (update_nn) {
    //     for (auto u : g.matrix[v]) {
    //         if (g.nodeset.fast_contain(u))
    //             ++num_neighbors_of[clq][u];
    //     }
    // }
}
// insert v into the col^th color. assumes it fits. Puts vertices
// added from candidates[i] into diff
template <class adjacency_struct>
void clique_finder<adjacency_struct>::insert_color(int v, int clq, bitset& diff)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    diff.copy(g.matrix[v]);
    diff.setminus_with(candidates[clq]);
    candidates[clq].union_with(g.matrix[v]);
}

template <class graph_struct>
void degeneracy_finder<graph_struct>::clear()
{
    order.clear();
}

// heuristically find a set of cliques and return the size of the
// largest
template <class graph_struct>
void degeneracy_finder<graph_struct>::degeneracy_ordering()
{
    order.reserve(g.size());
    // if (g.size() == g.capacity()) {
    for (auto v : g.nodes) {
        auto vd = g.matrix[v].size();
        if (vd >= buckets.size())
            buckets.resize(vd + 1);
        buckets[vd].push_front(v);
        degrees[v] = vd;
        iterators[v] = buckets[vd].begin();
        ordered[v] = false;
    }
    // } else {
    //     gc::bitset actual_neighbors(0, g.capacity() - 1, gc::bitset::empt);
    //     for (auto v : g.nodes) {
    //         actual_neighbors.copy(g.matrix[v]);
    //         actual_neighbors.intersect_with(g.nodeset);
    //         auto vd = actual_neighbors.size();
    //         if (vd >= buckets.size())
    //             buckets.resize(vd + 1);
    //         buckets[vd].push_front(v);
    //         degrees[v] = vd;
    //         iterators[v] = buckets[vd].begin();
    //         ordered[v] = false;
    //     }
    // }

    while (true) {
        size_t i{0};
        for (; i != buckets.size(); ++i)
            if (!buckets[i].empty())
                break;
        if (i == buckets.size())
            break;

        if (degeneracy < static_cast<int>(i)) {
            degeneracy = i;
            core.push_back(end(order));
            core_degree.push_back(degeneracy);
        }
        // d = std::max(d, static_cast<int>(i));

        auto v = buckets[i].back();
        order.push_back(v);
        buckets[i].pop_back();
        ordered[v] = true;

        for (auto u : g.matrix[v]) {
            if (!g.nodeset.fast_contain(u) or ordered[u])
                continue;
            auto& ud = degrees[u];
            buckets[ud].erase(iterators[u]);
            --ud;
            buckets[ud].push_front(u);
            iterators[u] = buckets[ud].begin();
        }
    }

    core.push_back(end(order));
    core_degree.push_back(g.size());
}

template <class graph_struct>
void degeneracy_finder<graph_struct>::display_ordering()
{
    std::cout << "Degeneracy ordering : " << std::endl;
    ;
    for (auto v : order)
        std::cout << v + 1 << "(" << degrees[v] << ")" << std::endl;
    std::cout << std::endl;
}

namespace detail
{
    template <class adjacency_struct> struct brelaz_state {
        clique_finder<adjacency_struct> cf;

        // degrees are valid only if the corresponding bit is not set in dirty
        mutable std::vector<int> degrees;
        mutable bitset dirty;

        // current set of uncolored nodes
        bitset nodes;

        // saturation
        std::vector<int> saturation;

        // temp
        mutable bitset util_set;

        void remove(int v)
        {
            nodes.fast_remove(v);
            dirty.union_with(cf.g.matrix[v]);
        }
        int degree(int v) const
        {
            if (dirty.fast_contain(v)) {
                util_set.copy(cf.g.matrix[v]);
                util_set.intersect_with(nodes);
                degrees[v] = util_set.size();
                dirty.fast_remove(v);
            }
            return degrees[v];
        }

        brelaz_state(const graph<adjacency_struct>& g)
            : cf(g)
            , degrees(cf.g.capacity(), 0)
            , dirty(0, cf.g.capacity(), bitset::full)
            , nodes(cf.g.nodeset)
            , saturation(cf.g.capacity(), 0)
            , util_set(0, cf.g.capacity(), bitset::empt)
        {
        }
    };

    template <class adjacency_struct> struct saturation_gt {
        const brelaz_state<adjacency_struct>& bs;

        bool operator()(int u, int v)
        {
            int satlt = bs.saturation[u] - bs.saturation[v];
            if (satlt > 0)
                return true;
            if (satlt < 0)
                return false;
            return bs.degree(u) > bs.degree(v);
        }
    };
} // namespace detail

template <class adjacency_struct>
std::vector<int> brelaz_color(
    const graph<adjacency_struct>& g, const bool randomized = false)
{
    detail::brelaz_state<adjacency_struct> state{g};
    auto& cf = state.cf;
    cf.clear();
    std::vector<int> solution(cf.g.capacity(), -1);
    Heap<detail::saturation_gt<adjacency_struct>> sheap(
        detail::saturation_gt<adjacency_struct>{state});

    if (randomized) {
        std::random_device rd;
        std::mt19937 s(rd());
        std::vector<int> order;
        for (auto v : cf.g.nodes)
            order.push_back(v);
        std::shuffle(begin(order), end(order), s);
        for (auto v : order)
            sheap.insert(v);

    } else
        for (auto v : cf.g.nodes)
            sheap.insert(v);

    // std::cout << std::endl << cf.num_cliques << std::endl;

    while (!sheap.empty()) {
        int v = sheap.removeMin();
        state.remove(v); // O(N/64) on bitsets, O(|N(v)|) on lists

        bool found{false};
        for (int i = 0; i != cf.num_cliques; ++i) {
            if (!cf.candidates[i].fast_contain(v)) {
                found = true;
                state.util_set.clear();

                // for (int j = 0; j != cf.num_cliques; ++j) {
                //      std::cout << "v[" << j << "] = " ;
                //      for( auto u : cf.candidates[j] )
                //              std::cout << " " << u;
                //      std::cout << std::endl;
                // }

                cf.insert_color(v, i, state.util_set);

                // std::cout << "color " << v << " with " << i << " (" <<
                // g.matrix[v].size() << "):"; for( auto u : g.matrix[v] ) {
                //      if( state.nodes.contain(u) ) {
                //              std::cout << " " << u;
                //      }
                // }
                // std::cout << std::endl ;
                // std::cout << std::endl ;

                state.util_set.intersect_with(state.nodes);
                for (auto u : state.util_set) {
                    ++state.saturation[u];
                    sheap.update(u);
                }
                solution[v] = i;
                break;
            }
        }
        if (!found) {

            cf.new_color();

            // for (int j = 0; j != cf.num_cliques; ++j) {
            //      std::cout << "v[" << j << "] = " ;
            //      for( auto u : cf.candidates[j] )
            //              std::cout << " " << u;
            //      std::cout << std::endl;
            // }

            cf.insert_color(v, cf.num_cliques - 1, state.util_set);
            solution[v] = cf.num_cliques - 1;

            // std::cout << "color " << v << " with " << solution[v] << " (" <<
            // g.matrix[v].size() << "):"; for( auto u : g.matrix[v] ) {
            //      if( state.nodes.contain(u) ) {
            //              std::cout << " " << u;
            //      }
            // }
            // std::cout << std::endl ;
            // std::cout << std::endl ;
        }
    }

    return solution;
}

// std::ostream& display_adjacency3(std::ostream& os,
// gc::basic_graph<gc::vertices_vec>& g)
//{
//      os << "Displaying adjacency\n";
//    for (auto it : g.nodes) {
//        os << "Vertex " << it+1 << "   Neighbors : ( ";
//              for (std::vector<int>::const_iterator itN =
//              g.matrix[it].vertices.begin(); itN <
//              g.matrix[it].vertices.end(); ++itN) {
//                      os << *itN+1 << " ";
//              }
//          os << ")";
//              os << " degree =" << g.matrix[it].vertices.size() << "\n";
//    }
//    return os;
//}

template <class adjacency_struct>
void minimal_triangulator<adjacency_struct>::section_node(int v)
{
    //      // PRINT
    //      std::cout << "Removing node :" << v+1 << std::endl;
    g_section.remove_node(v);
    // Remove o from the adjacency matrix of its neighbors
    for (auto u : g_section.matrix[v]) {
        std::vector<int>::iterator itr;
        itr = std::find(
            g_section.matrix[u].begin(), g_section.matrix[u].end(), v);
        if (itr != std::end(g_section.matrix[u])
            || (*g_section.matrix[u].end() == v)) {
            g_section.matrix[u].erase(itr);
        }
    }
}

// TODO: Genericity for dynamic sparse
template <class adjacency_struct>
void minimal_triangulator<adjacency_struct>::elimination_game(
    const std::vector<int>& ordering)
{
    int num_fill(0);
    for (auto o : ordering) {
        //              // PRINT
        //              std::cout << o+1 << std::endl;
        if (g_section.matrix[o].size() < 1) {
            // remove o !
            section_node(o);
            continue;
        } else {
            // Check every pair n1n2 of neighbors of o in g_section
            for (std::vector<int>::iterator n1 = g_section.matrix[o].begin();
                 n1 != (g_section.matrix[o].end() - 1); ++n1) {
                //                              // PRINT
                //                              std::cout << "   " << *n1+1 <<
                //                              std::endl;
                for (std::vector<int>::iterator n2 = n1 + 1;
                     n2 != g_section.matrix[o].end(); ++n2) {
                    //                                      // PRINT
                    //                                      std::cout << " " <<
                    //                                      *n2+1 << std::endl;
                    // If n1n2 are not neighors add to g_section and g_filled
                    std::vector<int>::iterator it;
                    it = std::find(g_section.matrix[*n1].begin(),
                        g_section.matrix[*n1].end(), *n2);
                    if (it == std::end(g_section.matrix[*n1])
                        && (*g_section.matrix[*n1].end() != *n2)) {
                        g_filled.add_edge(*n1, *n2);
                        g_section.add_edge(*n1, *n2);
                        // Need to sort adjacency matrix in case vertices_vec
                        // are used
                        std::sort(g_section.nodes.begin(),
                            g_section.nodes.end(),
                            [](int a, int b) { return a < b; });
                        std::pair<int, int> edge = std::make_pair(*n1, *n2);
                        fill_edges.push_back(edge);
                        std::cout << "Adding edge: (" << *n1 + 1 << ","
                                  << *n2 + 1 << ")" << std::endl;
                        ++num_fill;
                    }
                }
            }

            // Removing node o in section graph
            section_node(o);

            //                      // DISPLAYING GRAPHS AND FILL EDGES
            //            // Display g_section
            //            std::cout << "Eliminated graph :" << std::endl;
            //            for (auto its : g_section.nodes) {
            //              std::cout << "Vertex " << its+1 << "   Neighbors : (
            //              ";
            //                for (std::vector<int>::const_iterator itN =
            //                g_section.matrix[its].vertices.begin(); itN <
            //                g_section.matrix[its].vertices.end(); ++itN) {
            //                      std::cout << *itN+1 << " ";
            //                }
            //                std::cout << ")";
            //                std::cout << " degree =" <<
            //                g_section.matrix[its].vertices.size() <<
            //                std::endl;
            //            }

            //            // Display g_filled
            //            std::cout << "Filled graph :" << std::endl;
            //            for (auto itf : g_filled.nodes) {
            //              std::cout << "Vertex " << itf+1 << "   Neighbors : (
            //              ";
            //                for (std::vector<int>::const_iterator itN =
            //                g_filled.matrix[itf].vertices.begin(); itN <
            //                g_filled.matrix[itf].vertices.end(); ++itN) {
            //                      std::cout << *itN+1 << " ";
            //                }
            //                std::cout << ")";
            //                std::cout << " degree =" <<
            //                g_filled.matrix[itf].vertices.size() << std::endl;
            //                }
            //                      }
            //    }
            //    std::cout << "# fill edges: " << num_fill << std::endl;
            //    for (auto edge : fill_edges) {
            //        std::cout << "(" << edge.first+1 << "," << edge.second+1
            //        << ")        ";
            //    }
            //    std::cout << std::endl;
        }
    }
}

template <class adjacency_struct>
void minimal_triangulator<adjacency_struct>::stable_sort()
{
    for (auto u : g.nodes) {
        minimum_degree_ordering.push_back(u);
    }
    std::stable_sort(minimum_degree_ordering.begin(),
        minimum_degree_ordering.end(),
        [&](int a, int b) { return g.matrix[a].size() < g.matrix[b].size(); });
}

} // namespace gc

#endif
