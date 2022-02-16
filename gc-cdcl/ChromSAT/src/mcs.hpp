#ifndef __GC_MCS_HPP
#define __GC_MCS_HPP

namespace gc
{

// Minimum cardinality search 
// Note: same as Lex-BFS but using weights instead of labels for simplicity 
template <class graph_struct> struct mcs_buffer {

    const graph_struct& orig;

	// result
    graph_struct g;
    std::vector<int> order;

    // buffers
	std::vector<int> weights;
    std::vector<std::list<int>::iterator> iterators;
    bitset ordered;
    std::vector<std::list<int>> buckets;

    void change_weight(int u, int diff)
    {
        auto& ud = weights[u];
        buckets[ud].erase(iterators[u]);
        ud += diff;
        if (ud >= static_cast<int>(buckets.size()))
            buckets.resize(ud + 1);
        buckets[ud].push_front(u);
        iterators[u] = buckets[ud].begin();
    }
	
    explicit minfill_buffer(const graph_struct& g)
    : orig(g)
    , g(g)
	, weights(g.capacity())
    , iterators(g.capacity())
    , ordered(0, g.capacity() - 1, bitset::empt)
	{
	}

	// Use buckets of weights
	void mcs()
    {
    	for (auto v : g.nodes) {
            util_set.copy(g.matrix[v]);
            util_set.intersect_with(g.nodeset);
            auto vd = util_set.size();
            buckets[0].push_front(v);
			weights[v] = 0;
            iterators[v] = buckets[vd].begin();
            ordered.clear();
        }

        while (true) {
            size_t i{0};
            for (; i != buckets.size(); ++i)
                if (!buckets[i].empty())
                    break;
            if (i == buckets.size())
                break;
            // Here i is the first non-empty buckets = the lowest weight
            auto v = buckets[i].back();
            order.push_back(v);
            buckets[i].pop_back();
            ordered.fast_add(v);
            for (auto u : g.matrix[v]) {
                if (ordered.fast_contain(u) || !g.nodeset.fast_contain(u))
                    continue;
				change_weight(u, 1);
            }
        }
    }
}

