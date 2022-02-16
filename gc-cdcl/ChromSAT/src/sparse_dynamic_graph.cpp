
#include "sparse_dynamic_graph.hpp"


// #define _VERIFY_MCGRAPH

using std::swap;

using namespace gc;


// remove d from the bag 'd'
void coloring::remove(const int y, const int d) {

    int idy = rank[y];

    // std::cout << "remove " << y << " @ " << idy << " from bag " << d <<
    // std::endl;

    // swap y with first[satur[y].size-1], and increment first[satur[y].size-1]
    int idf = first[d];
    int f = order[idf];

    // std::cout << " -> swap with " << f << " @ " << idf << std::endl;

    rank[y] = idf;
    rank[f] = idy;

    order[idf] = y;
    order[idy] = f;

    ++first[d];

    // exit(1);
}

template <class ForwardIt, class Compare>
void max_elements(
    ForwardIt first, ForwardIt last, std::vector<ForwardIt>& maxes, Compare lt)
{
    if (first == last)
        return;

    // ForwardIt largest = first;
    maxes.push_back(first);
    ++first;
    for (; first != last; ++first) {
        if (!lt(*first, *(maxes.back()))) {
            if (lt(*(maxes.back()), *first))
                maxes.clear();
            maxes.push_back(first);
        }
    }
}

void coloring::brelaz_color(dyngraph& g, const int randomized = 0)
{
    if (g.nodes.empty())
        return;

    int iter = 0, c, d, x;

    rank.resize(g.capacity);
    color.resize(g.capacity, -1);
    satur.resize(g.capacity);

    //
    order.clear();

    // nodes are stored by non-increasing saturation degree
    for (auto v : g.nodes) {
        // rank[v] = order.size();
        order.push_back(v);
        // std::cout << v << " (" << g.matrix[v].size() << ")\n";
    }
    if (randomized == 1) {
        std::mt19937 s(rd());
        std::shuffle(begin(order), end(order), s);
        // std::cout << order[0] << std::endl;
    }

    x = 0;
    for (auto v : order)
        rank[v] = x++;

    first.clear();
    first.push_back(0); // every node has saturation degree 0

    std::vector<std::vector<int>::iterator> maxes;

    do {
        if (g.nodes.empty())
            break;

        // get the highest saturation degree
        d = satur[order[iter++]].size();

        // remove the unused pointers
        while (first.size() > d + 1)
            first.pop_back();

        // find out the vertex of highest degree among thos of
        // maximum saturation
        if (randomized == 2) {
            max_elements(begin(order) + first[d],
                (d > 0 ? begin(order) + first[d - 1] : end(order)), maxes,
                [&](const int x_, const int y_) {
                    return (g.degree(x_) < g.degree(y_));
                });
            x = *(maxes[rand() % maxes.size()]);
            maxes.clear();
        } else {
            x = *std::max_element(begin(order) + first[d],
                (d > 0 ? begin(order) + first[d - 1] : end(order)),
                [&](const int x_, const int y_) {
                    return (g.degree(x_) < g.degree(y_));
                });
        }
				// x = *(begin(order) + first[d]);

        // remove x from the partition of nodes with saturation
        // degree d
        remove(x, d);

        // use the first possible color for x
        c = satur[x].min();
        color[x] = c;

        // remove x from the graph
        g.rem_node(x);

        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (satur[y].remove(c)) {
                // new highest saturation degree, new pointer
                if (first.size() <= satur[y].size())
                    first.push_back(*rbegin(first));

                // move y one partition up in the saturation degree
                // list
                remove(y, satur[y].size() - 1);
            }

    } while (true);
    // std::cout << std::endl << std::endl;
}



dyngraph::dyngraph(const size_t n)
{
    capacity = n;

    nodes.reserve(capacity);
    nodes.fill();

    nodeset.initialise(0, capacity - 1, gc::bitset::full);

    matrix.resize(capacity);
    nb_index.resize(capacity);

    num_edges = 0;
}

dyngraph::dyngraph(const dyngraph& g)
{
    capacity = g.capacity;

    nodeset.initialise(0, capacity - 1, gc::bitset::empt);
    for (auto v : g.nodes) {
        declare_node(v);
    }

    matrix.resize(capacity);
    nb_index.resize(capacity);

    num_edges = g.num_edges;
    nodeset.copy(g.nodeset);
    for (unsigned i = 0; i < g.ranks.size(); ++i) {
        edges.push_back(g.edges[i]);
        ranks.push_back(g.ranks[i]);
    }

    for (int x = 0; x < capacity; ++x) {

        for (auto it = begin(g.matrix[x]); it != end(g.matrix[x]); ++it)
            matrix[x].push_back(*it);

        for (auto it = begin(g.nb_index[x]); it != end(g.nb_index[x]); ++it)
            nb_index[x].push_back(*it);
    }
}

size_t dyngraph::size() const { return nodes.size(); }

bool dyngraph::null() const { return nodes.size() == 0; }

bool dyngraph::empty() const { return num_edges == 0; }

bool dyngraph::full() const
{
    return nodes.size() * (nodes.size() - 1) == 2 * num_edges;
}

double dyngraph::get_density() const
{
    return (double)num_edges / (double)(size() * (size() - 1) / 2);
}

void dyngraph::sort(bool non_decreasing)
{
    // verify("before sort");

    std::vector<int> sorted;
    for (int i = 0; i < size(); ++i) {
        sorted.push_back(nodes[i]);
    }

    std::sort(sorted.begin(), sorted.end(),
        [&](const int x, const int y) { return (degree(x) > degree(y))^non_decreasing; });

    std::vector<int> srank(capacity);
    for (int i = 0; i < size(); ++i) {
        srank[sorted[i]] = i;
    }

    ranks.clear();
    num_edges = 0;
    std::vector<Edge> old_edges;
    std::swap(edges, old_edges);

    for (auto x = 0; x < capacity; ++x) {
        matrix[x].clear();
        nb_index[x].clear();
    }

    for (auto p : old_edges) {
        add_edge(srank[p[0]], srank[p[1]]);
    }

    // verify("after sort");
}


void dyngraph::clear()
{
    int x;
    while (!nodes.empty()) {
        x = nodes.back();
        nodes.pop_back();
        matrix[x].clear();
        nb_index[x].clear();
    }
    ranks.clear();
}

void dyngraph::declare_node(const int x)
{
    nodes.add(x);
    nodeset.add(x);

    // if(x >= capacity) {
    //              capacity = x+1;
    //                      matrix.resize(capacity);
    //                      nb_index.resize(capacity);
    // }
}

void dyngraph::add_node(const int x)
{
#ifdef _VERIFY_MCGRAPH
    verify("before add node");
#endif

    declare_node(x);

    int i = degree(x), y, e, pos;
    num_edges += i;

    while (i--) {
        y = matrix[x][i];
        e = nb_index[x][i];

        // we add x at the back of y's neighbor list
        pos = (e & 1);

        // we change the position of x in y's neighbor list
        ranks[e / 2][1 - pos] = matrix[y].size();

        // points to the edge from y's perspective
        nb_index[y].push_back(e ^ 1);

        // add x in y's neighbors
        matrix[y].push_back(x);

    }

#ifdef _VERIFY_MCGRAPH
    verify("after add node");
#endif
}

void dyngraph::rem_node(const int x)
{
    // std::cout << "remove " << x << std::endl;
#ifdef _VERIFY_MCGRAPH
    verify("before rem node");
#endif

    nodes.remove(x);
    nodeset.remove(x);
    auto i = degree(x);
    int y, z, rx, ex, posx, ey, posy;
    num_edges -= i;

    while (i--) {

        y = matrix[x][i];

        ex = nb_index[x][i];
        posx = ex & 1;

        // store the position of x in y's neighborhood
        rx = ranks[ex / 2][1 - posx];

        assert(edges[ex / 2][posx] == x);
        assert(edges[ex / 2][1 - posx] == y);
        assert(matrix[y][rx] == x);
        assert(ranks[ex / 2][posx] == i);

        // if(y == 73 or y == 74)
        // std::cout << " - from N(" << y << "): " << edges[ex / 2] << " where
        // it is at rank " << rx << "\n";

        // replace x by z
        z = matrix[y].back();
        matrix[y].pop_back();

        // if(y == 73 or y == 74)
        // std::cout << "replace edge " << edges[ex / 2] << " by edge " <<
        // edges[ nb_index[y].back() / 2 ] << std::endl;

        if (z != x) {
            matrix[y][rx] = z;

            // set the new position of z in y's neighborhood
            ey = nb_index[y].back();
            posy = (ey & 1);
            ranks[ey / 2][posy] = rx;

            //
            nb_index[y][rx] = ey;
        }

        nb_index[y].pop_back();
    }

#ifdef _VERIFY_MCGRAPH
    verify("after rem node");
#endif
}

bool dyngraph::has_node(int x) const { return nodeset.fast_contain(x); }

int dyngraph::add_edge(const int x, const int y)
{
    assert(nodes.contain(x));
    assert(nodes.contain(y));

    nb_index[x].push_back(2 * ranks.size());
    nb_index[y].push_back(2 * ranks.size() + 1);

    Edge r(matrix[x].size(), matrix[y].size());
    ranks.push_back(r);

    Edge e(x, y);
    edges.push_back(e);

    matrix[x].push_back(y);
    matrix[y].push_back(x);

    ++num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after add edge");
#endif

    return ranks.size() - 1;
}

void dyngraph::rem_edge(const int i)
{
    Edge edge = edges[i];
    rem_edge(edge[0], edge[1], i);
}
void dyngraph::rem_edge(const int x, const int y, const int e)
{

#ifdef _VERIFY_MCGRAPH
    assert(nodes.contain(x));
    assert(nodes.contain(y));
#endif

    int ry = ranks[e][0];
    int rx = ranks[e][1];

    int sx = matrix[y].back();
    matrix[y].pop_back();

    int ey = nb_index[y].back();
    nb_index[y].pop_back();
    if (x != sx) {
        int posy = (ey & 1);
        ranks[ey / 2][posy] = rx;

        matrix[y][rx] = sx;
        nb_index[y][rx] = ey;
    }

    int sy = matrix[x].back();
    matrix[x].pop_back();

    int ex = nb_index[x].back();
    nb_index[x].pop_back();
    if (y != sy) {
        int posx = (ex & 1);

        ranks[ex / 2][posx] = ry;

        matrix[x][ry] = sy;
        nb_index[x][ry] = ex;
    }

    Edge ez = edges.back();
    edges.pop_back();
    Edge rz = ranks.back();
    ranks.pop_back();

    if (static_cast<unsigned>(e) != edges.size()) {

        int s;
        for (int i = 0; i < 2; ++i) {
            s = (nb_index[ez[i]][rz[i]] & 1);
            nb_index[ez[i]][rz[i]] = 2 * e + s;
        }

        edges[e] = ez;
        ranks[e] = rz;
    }

    --num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after rem edge");
#endif
}

void dyngraph::maximal_matching(
    std::vector<int>& matching, int& nmatch, std::vector<int>& ranklist)
{
    ranklist.clear();
    for (int i = 0; i < size(); i++)
        ranklist.push_back(i);

    // std::random_device rd;
    // std::mt19937 g(rd());
    //     std::shuffle(ranklist.begin(), ranklist.end());

    int u, v;
    matching.assign(capacity, -1);
    nmatch = 0;
    for (auto i = 0; i < size(); i++) {
        u = nodes[ranklist[i]];
        if (matching[u] == -1) {
            for (size_t j = 0; j < matrix[u].size(); j++) {
                v = matrix[u][j];
                if (matching[v] == -1) {
                    matching[u] = v;
                    matching[v] = u;
                    nmatch++;
                    break;
                }
            }
        }
    }
}

std::ostream& dyngraph::display(std::ostream& os) const
{

    for (auto i = 0; i < size(); ++i) {
        int x = nodes[i];

        os << x << ": ";
        vecdisplay(matrix[x], os);
        os << " (" << degree(x) << ")" << std::endl;

        // os << "   [";
        // for(auto e : nb_index[x])
        //              os << " " << (e/2);
        // os << " ]\n";
    }
    // for(auto e : edges)
    //              std::cout << e << std::endl;

    return os;
}

void dyngraph::print_dimacs(std::ostream& os) const
{
    os << "c Generated by minicsp" << std::endl
       << "p edge " << size() << " " << edges.size() << std::endl;

    // not true during search
    assert(static_cast<int>(edges.size()) == num_edges);

    for (auto edge : edges) {
        os << "e " << (edge[0] + 1) << " " << (edge[1] + 1) << std::endl;
    }
}

void dyngraph::verify(const char* msg)
{
    std::cout << "check dyngraph consistency\n";

    for (unsigned i = 0; i < edges.size(); ++i) {
        Edge e = edges[i];
        int x = e[0];
        int y = e[1];

        Edge r = ranks[i];
        int ry = r[0];
        int rx = r[1];

        // std::cout << "rank of " << y << " N(" << x << ") = " << ry <<
        // std::endl; for()

        if (nodes.contain(y) && matrix[x][ry] != y) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(matrix[x], std::cout);
            std::cout << " is not " << y << std::endl;
            assert(0);
        }
        if (nodes.contain(y) && static_cast<unsigned>(nb_index[x][ry]) != 2 * i) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(nb_index[x], std::cout);
            std::cout << " is not " << (2 * i) << std::endl;
            assert(0);
        }

        if (nodes.contain(x) && matrix[y][rx] != x) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << rx << "-th element of ";
            vecdisplay(matrix[y], std::cout);
            std::cout << " is not " << x << std::endl;
            assert(0);
        }
        if (nodes.contain(x)
            && static_cast<unsigned>(nb_index[y][rx]) != 2 * i + 1) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(nb_index[y], std::cout);
            std::cout << " is not " << (2 * i + 1) << std::endl;
            assert(0);
        }
    }

    int ecount = 0;
    for (int i = 0; i < size(); ++i) {
        int x = nodes[i];

        for (unsigned j = 0; j < matrix[x].size(); ++j) {

            ++ecount;
            int y = matrix[x][j];
            int ey = nb_index[x][j];

            Edge e = edges[ey / 2];

            int posx = (ey & 1);

            if(e[posx] != x) {
                std::cout << msg << " " << j << "-th neighbor of " << x
                          << " is " << y << " but " << x << " is not the "
                          << posx << "-th node of edge " << e << std::endl;
                assert(0);
            }
            if(e[1 - posx] != y) {
                std::cout << msg << " " << j << "-th neighbor of " << x
                          << " is " << y << " but " << y << " is not the "
                          << 1 - posx << "-th node of edge " << e << std::endl;
                assert(0);
            }

            Edge r = ranks[ey / 2];

            if(static_cast<unsigned>(r[posx]) != j) {
                std::cout << msg << " " << x << "'s rank in N(" << x << ") is "
                          << r[posx] << ", should be " << j << std::endl;
                vecdisplay(nb_index[y], std::cout);
                std::cout << std::endl;
                assert(0);
            }
        }
    }
    ecount /= 2;

    if (ecount != num_edges) {
        std::cout << " " << msg << ": num_edges = " << num_edges
                  << " real count = " << ecount << std::endl
            ;

        assert(0);
    }
}

// void dyngraph::brelaz_color(coloring& col) {
//              if(node.empty()) return;
//
//              // first order the nodes by degree
//              col.degree.resize(node.size());
//              for(auto v : node) {
//                              col.order.push_back(v);
//                              col.degree[v] = degree(v);
//              }
//
//                  std::sort(col.order.begin(), col.order.end(),
//                      [&](const int x, const int y) { return (degree(x) >
//                      degree(y)); });
//
//
//              for(auto v : col.order) {
//                      std::cout << v << " " << degree(v) << std::endl;
//              }
//
//              int c;
//              int x = order[0]; // the next vertex to color
//              int num_colors = 0; // the number of colors used so far
//
//              col.satur.resize(node.size());
//              col.first.push_back(0); // every node has saturation degree 0
//
//              do {
//                              // forbidden_col.clear();
//                              // for(int i=neighbor[x].size();
//                              i<col.degree[i]; ++i) {
//                              // forbidden_col.add(col.color[neighbor[x][i]]);
//                              // }
//                              // for(c=0; c<num_colors; ++c)
//                              //              if(!forbidden_col.contain(c))
//                              break; c = col.satur[x].get(); col.color[x] = c;
//                              for(auto y : neighbor[x])
//                                              if(col.satur[y].add(c)) {
//                                                      if(col.first.size() <=
//                                                      col.satur[y].size)
//                                                      // swap
//
//                                              }
//
//
//              } while( !node.empty() );
//
//
//
// }

std::ostream& gc::operator<<(std::ostream& os, const dyngraph& x)
{
    return x.display(os);
}

std::ostream& gc::operator<<(std::ostream& os, const dyngraph* x)
{
    return (x ? x->display(os) : os);
}
