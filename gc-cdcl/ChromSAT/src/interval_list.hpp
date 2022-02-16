

#include <assert.h>
#include <iomanip>
#include <vector>




#ifndef __INTERVAL_LIST_HPP
#define __INTERVAL_LIST_HPP





namespace gc
{

class interval_list
{
	/*
	- next implements a linked list of indices for inf/sup the end of the list is signaled by "-1"
	- freed stores unused indices in order to optimise space
	- inf/sup store the lower and upper bounds of the interval, respectively
	
	*/
	

public:
    static const int BIG = 0x1000000;

    std::vector<int> inf;
    std::vector<int> sup;
    std::vector<int> next;
    std::vector<int> freed;

    size_t size_;

    inline size_t size() { return size_; }

    interval_list();

    bool contain(const int x);
    bool remove(const int x);
    int min() const;
    void fill();

    void initialise(const int ub) {} // dummy

    void check_consistency();

    std::ostream& display(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const interval_list& x);

std::ostream& operator<<(std::ostream& os, const interval_list* x);

} // namespace gc

#endif // __INTERVAL_LIST_HPP
