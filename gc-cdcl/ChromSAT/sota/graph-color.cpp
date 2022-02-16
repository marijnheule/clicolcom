/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Contributing authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Christian Schulte, 2004
 *     Stefano Gualandi, 2013
 *
 *  Last modified:
 *     $Date$ by $Author$
 *     $Revision$
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <gecode/driver.hh>
#include <gecode/int.hh>

#include <intstack.hpp>
#include <bitset.hpp>
#include <graph.hpp>
#include <dimacs.hpp>

using namespace Gecode;


/**
 * \name %Graph specification for graph coloring
 *
 * The edges are described by an array of integers with even number
 * of elements, terminated by the elements -1,-1.
 * The cliques are described by an array of integers, where the first
 * integer gives the size of the clique, the following elements are
 * nodes for each clique. The cliques are terminated by -1 for clique size
 *
 * \relates GraphColor
 */
//@{
/// %Graph specification
class GraphColorSpec
{
public:
    const int n_v; ///< Number of nodes
    std::vector<int> e; ///< Edges
    std::vector<int> c; ///< Cliques
    GraphColorSpec(const int n_v0, std::vector<int> e0, std::vector<int> c0)
        : n_v(n_v0)
        , e(e0)
        , c(c0)
    {
    }
};

GraphColorSpec read_dimacs(const std::string& fname)
{
    gc::graph g;
    std::vector<int> edges;
    dimacs::read_graph(fname.c_str(), [&](int nv, int) { g = gc::graph{nv}; },
        [&](int u, int v) {
            if (u != v) {
                g.add_edge(u - 1, v - 1);
                edges.push_back(u - 1);
                edges.push_back(v - 1);
            }
        },
        [&](int, gc::weight) {});
    edges.push_back(-1);
    edges.push_back(-1);
    std::cout << "Read graph" << std::endl;
    gc::clique_finder cf(g);
    cf.find_cliques(g.nodes);
    std::vector<int> cliques;
    std::cout << "found " << cf.clique_sz.size() << " cliques\n";
    for (auto& clq : cf.cliques) {
      cliques.push_back(clq.size());
      for (auto v : clq) {
        cliques.push_back(v);
      }
    }
    cliques.push_back(-1);
    GraphColorSpec s(g.capacity(), edges, cliques);
    return s;
}

/**
 * \brief %Example: Clique-based graph coloring
 *
 * \ingroup Example
 *
 */
class GraphColor : public IntMinimizeScript
{
private:
    const GraphColorSpec& g;
    /// Color of nodes
    IntVarArray v;
    /// Number of colors
    IntVar m;

public:
    /// Model variants
    enum {
        MODEL_NONE, ///< No lower bound
        MODEL_CLIQUE ///< Use maximal clique size as lower bound
    };
    /// Branching to use for model
    enum {
        BRANCH_DEGREE, ///< Choose variable with largest degree
        BRANCH_SIZE, ///< Choose variablee with smallest size
        BRANCH_BRELAZ, ///< Choose variablee with smallest size
        BRANCH_DEGREE_SIZE, ///< Choose variable with largest degree/size
        BRANCH_AFC_SIZE, ///< Choose variable with largest afc/size
        BRANCH_ACTION_SIZE ///< Choose variable with smallest size/action
    };
    /// Symmetry variants
    enum {
        SYMMETRY_NONE, ///< Simple symmetry
        SYMMETRY_LDSB ///< Use LDSB for value symmetry breaking
    };
    /// The actual model
    GraphColor(const InstanceOptions& opt)
        : IntMinimizeScript(opt)
        , g(read_dimacs(opt.instance()))
        , v(*this, g.n_v, 0, g.n_v - 1)
        , m(*this, 0, g.n_v - 1)
    {
        rel(*this, v, IRT_LQ, m);
        for (int i = 0; g.e[i] != -1; i += 2)
            rel(*this, v[g.e[i]], IRT_NQ, v[g.e[i + 1]]);

        const int* c = g.c.data();
        for (int i = *c++; i--; c++)
            rel(*this, v[*c], IRT_EQ, i);
        while (*c != -1) {
            int n = *c;
            IntVarArgs x(n);
            c++;
            for (int i = n; i--; c++)
                x[i] = v[*c];
            distinct(*this, x, opt.ipl());
            if (opt.model() == MODEL_CLIQUE)
                rel(*this, m, IRT_GQ, n - 1);
        }
        // Branching on the number of colors
        branch(*this, m, INT_VAL_MIN());
        if (opt.symmetry() == SYMMETRY_NONE) {
            // Branching without symmetry breaking
            switch (opt.branching()) {
            case BRANCH_SIZE:
                branch(*this, v, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
                break;
            case BRANCH_BRELAZ:
                branch(*this, v,
                    tiebreak(INT_VAR_SIZE_MIN(), INT_VAR_DEGREE_MAX()),
                    INT_VAL_MIN());
                break;
            case BRANCH_DEGREE:
                branch(*this, v,
                    tiebreak(INT_VAR_DEGREE_MAX(), INT_VAR_SIZE_MIN()),
                    INT_VAL_MIN());
                break;
            case BRANCH_DEGREE_SIZE:
                branch(*this, v, INT_VAR_DEGREE_SIZE_MAX(), INT_VAL_MIN());
                break;
            case BRANCH_AFC_SIZE:
                branch(
                    *this, v, INT_VAR_AFC_SIZE_MAX(opt.decay()), INT_VAL_MIN());
                break;
            case BRANCH_ACTION_SIZE:
                branch(*this, v, INT_VAR_ACTION_SIZE_MAX(opt.decay()),
                    INT_VAL_MIN());
                break;
            }
        } else { // opt.symmetry() == SYMMETRY_LDSB
            // Branching while considering value symmetry breaking
            // (every permutation of color values gives equivalent
            // solutions)
            Symmetries syms;
            syms << ValueSymmetry(IntArgs::create(g.n_v, 0));
            switch (opt.branching()) {
            case BRANCH_SIZE:
                branch(*this, v, INT_VAR_SIZE_MIN(), INT_VAL_MIN(), syms);
                break;
            case BRANCH_BRELAZ:
                branch(*this, v,
                    tiebreak(INT_VAR_SIZE_MIN(), INT_VAR_DEGREE_MAX()),
                    INT_VAL_MIN(), syms);
                break;
            case BRANCH_DEGREE:
                branch(*this, v,
                    tiebreak(INT_VAR_DEGREE_MAX(), INT_VAR_SIZE_MIN()),
                    INT_VAL_MIN(), syms);
                break;
            case BRANCH_DEGREE_SIZE:
                branch(
                    *this, v, INT_VAR_DEGREE_SIZE_MAX(), INT_VAL_MIN(), syms);
                break;
            case BRANCH_AFC_SIZE:
                branch(*this, v, INT_VAR_AFC_SIZE_MAX(opt.decay()),
                    INT_VAL_MIN(), syms);
                break;
            case BRANCH_ACTION_SIZE:
                branch(*this, v, INT_VAR_ACTION_SIZE_MAX(opt.decay()),
                    INT_VAL_MIN(), syms);
                break;
            }
        }
    }
    /// Cost function
    virtual IntVar cost(void) const { return m; }
    /// Constructor for cloning \a s
    GraphColor(GraphColor& s)
        : IntMinimizeScript(s)
        , g(s.g)
    {
        v.update(*this, s.v);
        m.update(*this, s.m);
    }
    /// Copying during cloning
    virtual Space* copy(void) { return new GraphColor(*this); }
    /// Print the solution
    virtual void print(std::ostream& os) const
    {
        int c = m.val() + 1;
        os << "\tm = " << c << std::endl << "\tv[] = {";
        for (int i = 0; i < v.size(); i++) {
            os << v[i] << ", ";
            if ((i + 1) % 15 == 0)
                os << std::endl << "\t       ";
        }
        os << "};" << std::endl;
    }
};

/** \brief Main-function
 *  \relates GraphColor
 */
int main(int argc, char* argv[])
{
    InstanceOptions opt("GraphColor");
    opt.ipl(IPL_DOM);
    opt.solutions(0);

    opt.model(GraphColor::MODEL_NONE);
    opt.model(GraphColor::MODEL_NONE, "none", "no lower bound");
    opt.model(GraphColor::MODEL_CLIQUE, "clique",
        "use maximal clique size as lower bound");

    opt.branching(GraphColor::BRANCH_DEGREE);
    opt.branching(GraphColor::BRANCH_DEGREE, "degree");
    opt.branching(GraphColor::BRANCH_SIZE, "size");
    opt.branching(GraphColor::BRANCH_BRELAZ, "brelaz");
    opt.branching(GraphColor::BRANCH_DEGREE_SIZE, "degree-size");
    opt.branching(GraphColor::BRANCH_AFC_SIZE, "afc-size");
    opt.branching(GraphColor::BRANCH_ACTION_SIZE, "action-size");

    opt.symmetry(GraphColor::SYMMETRY_NONE);
    opt.symmetry(GraphColor::SYMMETRY_NONE, "none");
    opt.symmetry(
        GraphColor::SYMMETRY_LDSB, "ldsb", "use value symmetry breaking");

    opt.parse(argc, argv);
    Script::run<GraphColor, BAB, InstanceOptions>(opt);
    return 0;
}

// STATISTICS: example-any
