#ifndef __GC_VARMAP_HPP
#define __GC_VARMAP_HPP

#include "minicsp/core/solvertypes.hpp"
#include <sparsehash/dense_hash_map>


namespace gc
{

struct varmap {
    using inner_map = google::dense_hash_map<int, int>;
    using hmap = google::dense_hash_map<int, inner_map>;
    hmap vars;

    struct proxy {
        const inner_map* m;
        minicsp::Var operator[](int u) const {
            if (!m)
                return minicsp::var_Undef;
            auto j = m->find(u);
            if (j == m->end())
                return minicsp::var_Undef;
            return j->second;
        }
    };

    bool contain(int v, int u) const
    {
        auto i = vars.find(v);
        if (i == vars.end()) {
            return false;
        }
        return i->second.find(u) != i->second.end();
    }

    proxy operator[](int v) const {
        auto i = vars.find(v);
        if (i == vars.end()) {
            return proxy{nullptr};
        }
        return proxy{&(i->second)};
    }

    varmap() { vars.set_empty_key(-1); }

    template<typename I, typename S>
    explicit varmap(I begin, S end)
    {
        vars.set_empty_key(-1);
        for (I i = begin; i != end; ++i)
            vars[*i].set_empty_key(-1);
    }

    varmap(const varmap&) = default;
    varmap(varmap&&) = default;
    varmap& operator=(const varmap&) = default;
    varmap& operator=(varmap&&) = default;
};

}

#endif
