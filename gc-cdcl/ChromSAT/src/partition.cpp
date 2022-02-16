
#include "partition.hpp"

    /*!@name Constructors*/
    //@{
    partition::partition() {}

    int partition::size() { return bag.size(); }

    void partition::clear() { bag.clear(); }
    void partition::resize(const size_t n, const size_t m)
    {
        bag.resize(m);
        index_.resize(n, NOTIN);
    }

    const std::vector<int>& partition::operator[](const int i) const
    {
        return bag[i];
    }

    /*!@name List Manipulation*/
    //@{
    void partition::move(const int elt, const int from, const int to) {

        // std::cout << "move " << elt << " from " << from << " to " << to <<
        // "\n";

        bag[from][index_[elt]] = bag[from].back();
        index_[bag[from].back()] = index_[elt];
        bag[from].pop_back();
        add_elt(elt, to);
    }
    void partition::add_elt(const int elt, const int to)
    {

        // std::cout << "add " << elt  << " to (" << to << ")\n";
        // for(auto e : bag[to]) {
        //         // 	std::cout << " " << e;
        // 	assert(e != elt);
        // }
        //         // std::cout << std::endl;

        index_[elt] = bag[to].size();
        bag[to].push_back(elt);
    }
    bool partition::contain(const int elt, const int c)
    {
        auto i{index_[elt]};
        return (i < bag[c].size() and bag[c][i] == elt);
    }
    void partition::swap(const int a, const int b)
    {
        std::swap(bag[a], bag[b]);
    }
    void partition::remove(const int a)
    {
        swap(a, static_cast<int>(bag.size()) - 1);
        bag.pop_back();
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& partition::display(std::ostream& os) const
    {
        for (auto it{begin(bag)}; it != end(bag); ++it) {
            os << (it - begin(bag)) << ":";
            for (auto jt{begin(*it)}; jt != end(*it); ++jt) {
                os << " " << *jt;
            }
            os << std::endl;
        }
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const partition& x)
    {
        return x.display(os);
    }
