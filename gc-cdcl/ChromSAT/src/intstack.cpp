
#include <assert.h>

#include "intstack.hpp"

intstack::intstack(const size_t n)
{
    size_ = 0;
    reserve(n);
}

void intstack::reserve(const size_t n)
{
    while (list_.size() < n) {
        index_.push_back(list_.size());
        list_.push_back(list_.size());
    }
}

// void intstack::resize(const size_t n)
// {
//     while (list_.size() < n) {
//         index_.push_back(list_.size());
//         list_.push_back(list_.size());
//     }
// 		// while (list_.size() > n) {
// 		// 	    auto last = list_.back();
// 		// 	list_.pop_back();
// 		// 	auto idmax = index_[list_.size()];
// 		// 	list_[idmax] = last;
// 		// 	index_.pop_back();
// 		// }
// }

void intstack::save(size_t& stamp) { stamp = size_; }
void intstack::restore(const size_t stamp) { size_ = stamp; }

//@}

/*!@name Accessors*/
//@{

bool intstack::safe_contain(const int elt) const
{
    if (elt >= 0 && (size_t)elt < index_.size())
        return contain(elt);
    return false;
}

bool intstack::contain(const int elt) const { return index_[elt] < size_; }

size_t intstack::size() const { return size_; }

bool intstack::empty() const { return size_ == 0; }

int intstack::next(const int elt) const
{
    size_t idx = index_[elt] + 1;
    return (idx < size_ ? list_[idx] : elt);
}
int intstack::prev(const int elt) const
{
    size_t idx = index_[elt];
    return (idx > 0 ? list_[idx - 1] : elt);
}

int intstack::operator[](const size_t idx) const { return list_[idx]; }

int& intstack::operator[](const size_t idx) { return list_[idx]; }
//@}

/*!@name List Manipulation*/
//@{
std::vector<int>::iterator intstack::begin() { return list_.begin(); }
std::vector<int>::reverse_iterator intstack::rbegin()
{
    return list_.rend() + size_;
}

std::vector<int>::iterator intstack::end() { return list_.begin() + size_; }
std::vector<int>::reverse_iterator intstack::rend() { return list_.rend(); }

std::vector<int>::const_iterator intstack::begin() const
{
    return list_.begin();
}
std::vector<int>::const_reverse_iterator intstack::rbegin() const
{
    return list_.rend() + size_;
}

std::vector<int>::const_iterator intstack::end() const
{
    return list_.begin() + size_;
}
std::vector<int>::const_reverse_iterator intstack::rend() const
{
    return list_.rend();
}

std::vector<int>::iterator intstack::begin_not_in() { return end(); }
std::vector<int>::reverse_iterator intstack::rbegin_not_in()
{
    return list_.rend();
}

std::vector<int>::iterator intstack::end_not_in() { return list_.end(); }
std::vector<int>::reverse_iterator intstack::rend_not_in() { return rbegin(); }

std::vector<int>::const_iterator intstack::begin_not_in() const
{
    return end();
}
std::vector<int>::const_reverse_iterator intstack::rbegin_not_in() const
{
    return list_.rend();
}

std::vector<int>::const_iterator intstack::end_not_in() const
{
    return list_.end();
}
std::vector<int>::const_reverse_iterator intstack::rend_not_in() const
{
    return rend();
}

void intstack::fill() { size_ = list_.size(); }

void intstack::clear() { size_ = 0; }

// void intstack::set_size(const int s) { size_ = s; }

void intstack::safe_remove(const int elt)
{
    if (elt >= 0) {
        if (static_cast<size_t>(elt) >= list_.size()) {
            reserve(elt + 1);
        }
        remove(elt);
    }
}

void intstack::remove(const int elt)
{
    if (index_[elt] < size_)
        pull(elt);
}

void intstack::pull(const int elt)
{
    auto last = list_[--size_];
    index_[last] = index_[elt];
    list_[index_[elt]] = last;
    list_[size_] = elt;
    index_[elt] = size_;
}

void intstack::move_up(const int elt, const int idx_to)
{
    auto idx_from = index_[elt];

    assert(index_[elt] <= static_cast<size_t>(idx_to));

    auto last = list_[idx_to];
    index_[last] = idx_from;
    list_[idx_from] = last;
    list_[idx_to] = elt;
    index_[elt] = idx_to;
}

void intstack::pop_back() { --size_; }
void intstack::pop_head() { remove(list_[0]); }

int intstack::head() const { return list_[0]; }
int intstack::back() const { return list_[size_ - 1]; }

void intstack::safe_add(const int elt)
{
    if (elt >= 0) {
        if (static_cast<size_t>(elt) >= list_.size()) {
            reserve(elt + 1);
        }
        add(elt);
    }
}

void intstack::add(const int elt)
{
    if (index_[elt] >= size_)
        push(elt);
}

void intstack::push(const int elt)
{
    auto next = list_[size_];
    index_[next] = index_[elt];
    list_[index_[elt]] = next;
    index_[elt] = size_;
    list_[size_++] = elt;
}

int intstack::index(const int elt) const
{
		return index_[elt];
}
//@}

std::ostream& intstack::display(std::ostream& os) const
{
    os << "(";
    for (auto it = begin(); it < end(); ++it) {
        os << " " << *it;
    }
    os << " )";
    return os;
}

std::ostream& operator<<(std::ostream& os, const intstack& x)
{
    return x.display(os);
}

std::ostream& operator<<(std::ostream& os, const intstack* x)
{
    return (x ? x->display(os) : os);
}
