#include <iostream>

#include <algorithm>

#include "graph.hpp"

#ifndef __CLIQUESAMPLER_HPP
#define __CLIQUESAMPLER_HPP

namespace gc
{

static unsigned long x = 123456789, y = 362436069, z = 521288629;

unsigned long xorshf96(void)
{ // period 2^96-1
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

struct clique_sampler {

    std::vector<int> start_set;
    std::vector<int> cand_set;
    std::vector<int> probed;
    std::vector<std::vector<int>> domain;
    std::vector<int> buffer;

    gc::bitset nodeset;

    std::vector<int> clique;
    int lb;
    size_t probewidth;

    clique_sampler() {}

    template <class viterator>
    clique_sampler(viterator first, viterator last, const int n)
    {
        set_domain(first, last, n, true);
    }

    template <class viterator>
    void set_domain(
        viterator first, viterator last, const int n, const bool full)
    {
        nodeset.reinitialise(0, n - 1, 0);
        start_set.clear();

        if (full)
            nodeset.fill();

        for (auto vp{first}; vp != last; ++vp) {
            start_set.push_back(*vp);
            nodeset.fast_add(*vp);

#ifdef _DEBUG_SAMPLE
            std::cout << " " << (*vp);
#endif
        }

#ifdef _DEBUG_SAMPLE
        std::cout << std::endl;
#endif
    }

    template <class graph_struct, class viterator>
    int find_clique(graph_struct& g, const int l, viterator first,
        viterator last, const size_t basewidth, const size_t pw = 1)
    {
        lb = l;
        probewidth = pw;

        domain.resize(probewidth);

				// int count{0};
        for (auto vp{first}; vp != last; ++vp) {
            start_set.push_back(*vp);

						// if(count >= g.size()) {
						// 	if(g.nodes.contain(*vp)) {
						// 		std::cout << "GRR\n";
						// 		exit(1);
						// 	}
						// }
						// ++count;

#ifdef _DEBUG_SAMPLE
            std::cout << " " << (*vp);
#endif
        }

#ifdef _DEBUG_SAMPLE
        if (first != last)
            std::cout << std::endl;
#endif

        auto limit{std::min(basewidth, start_set.size())};

        for (auto i{0}; i < limit; ++i) {
            std::swap(start_set[i],
                start_set[i + (xorshf96() % (start_set.size() - i))]);
        }

        start_set.resize(limit);

        for (auto v : start_set) {
            
						clique.clear();
            clique.push_back(v);
						
            cand_set.clear();
            for (auto u : g.matrix[v])
                if (nodeset.fast_contain(u))
                    cand_set.push_back(u);
						probed = cand_set;

            auto p{probe(g)};
            // std::cout << " " << p ;
            // std::cout.flush();
            lb = std::max(lb, p);
        }
				// std::cout << std::endl;

        return lb;
    }

    template <class graph_struct> int probe(graph_struct& g)
    {

#ifdef _DEBUG_SAMPLE
        std::cout << "probe (" << lb
                  << "): " << print_container<std::vector<int>>(clique) << " + "
                  << print_container<std::vector<int>>(cand_set) << std::endl;
#endif

        if (clique.size() + cand_set.size() <= lb or cand_set.size() == 0)
            return clique.size();

        auto w{std::min(probewidth, cand_set.size())};

        int largest{0};
        // select 'probewidth' samples
        for (int i = 0; i < w; ++i) {
            std::swap(
                probed[i], probed[i + (xorshf96() % (probed.size() - i))]);

            auto v{cand_set[i]};
            domain[i].clear();
            for (auto u : g.matrix[v])
                if (nodeset.fast_contain(u))
                    domain[i].push_back(u);

            // Put intersection in the buffer
            buffer.clear();
            std::set_intersection(domain[i].begin(), domain[i].end(),
                cand_set.begin(), cand_set.end(), std::back_inserter(buffer));

            // Copy the buffer in domain[i]
            std::swap(buffer, domain[i]);

#ifdef _DEBUG_SAMPLE
            std::cout << v << " -> " << domain[i].size() << std::endl;
#endif

            if (domain[largest].size() < domain[i].size())
                largest = i;
        }

        clique.push_back(cand_set[largest]);

        std::swap(cand_set, domain[largest]);
        probed = cand_set;

        return probe(g);
    }
};

} // namespace gc

#endif // __CLIQUESAMPLER_HPP
