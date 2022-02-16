#include "vertices_vec.hpp"
#include "bitset.hpp"
#include "utils.hpp"
#include "minicsp/mtl/Alg.h"

// TODO
/* Copying : clear vector before ? use  = , std::copy, or else ?

*/

namespace gc
{

using bitset = BitSet;

/********************************************************/
//                                      VERTICES_VEC //
/********************************************************/
vertices_vec::vertices_vec(){}
vertices_vec::vertices_vec(vertices_vec const& v):vertices(v.vertices){}
vertices_vec::vertices_vec(const int lb, const int ub, const int p) { initialise(lb,ub,p); }
vertices_vec::~vertices_vec(){}


void vertices_vec::initialise(const int lb, const int ub, const int p) {
    if (p != 0) {
        for (int v = lb; v <= ub; ++v) {
            vertices.push_back(v);
        }
    }
}

void vertices_vec::copy(BitSet const& elts) {
    vertices.clear();
    for (auto v : elts)
        vertices.push_back(v);
}

void vertices_vec::canonize() {
    sort_by_number();
    vertices.erase(std::unique(begin(), end()), end());
}

std::vector<int>::const_iterator vertices_vec::find(const int elt) const
{
    auto i = std::lower_bound(begin(), end(), elt);
    if (i == end() or *i != elt)
        return end();
    return i;
}

std::vector<int>::iterator vertices_vec::find(const int elt)
{
    auto i = std::lower_bound(begin(), end(), elt);
    if (i == end() or *i != elt)
        return end();
    return i;
}

bool vertices_vec::fast_contain(const int elt) const {
    return find(elt) != end();
}

// Intersect with (Assuming both vectors are sorted)
void vertices_vec::intersect_with(vertices_vec const& v)
{
    // Put intersection in the buffer
    buffer.clear();
    std::set_intersection(vertices.begin(), vertices.end(), v.vertices.begin(),
        v.vertices.end(), std::back_inserter(buffer));

    // Copy the buffer in vertices
    // vertices.clear();
    swap(vertices, buffer);
}

void vertices_vec::intersect_with(bitset const& bs)
{
    using std::begin;
    using std::end;
    erase_if(vertices, [&](int elt) { return !bs.fast_contain(elt); });
    std::sort(vertices.begin(), vertices.end());
}

void vertices_vec::safe_intersect_with(vertices_vec & v)
{
    if (!v.is_sorted())
        v.sort_by_number();
    if (!is_sorted())
        sort_by_number();
    intersect_with(v);
}

// Union with another vertices_vec
void vertices_vec::union_with(vertices_vec const& v)
{
    // Put union in the buffer
    buffer.clear();
    std::set_union(vertices.begin(), vertices.end(), v.vertices.begin(),
        v.vertices.end(), std::back_inserter(buffer));

    // Copy the buffer in vertices
    // vertices.clear();
    vertices = buffer;
}

void vertices_vec::safe_union_with(vertices_vec & v)
{
    if (!v.is_sorted())
        v.sort_by_number();
    if (!is_sorted())
        sort_by_number();
    union_with(v);
}

// Union with a single vertex
void vertices_vec::union_with(int i)
{
    if (!is_sorted())
        sort_by_number();
    vertices.insert(std::upper_bound(vertices.begin(), vertices.end(), i), i);
}

// Difference with
void vertices_vec::difference_with(vertices_vec const& v)
{
    // Put difference in the buffer
    buffer.clear();
    std::set_difference(vertices.begin(), vertices.end(), v.vertices.begin(),
        v.vertices.end(), std::inserter(buffer, buffer.begin()));
    // Copy the buffer in vertices
    // vertices.clear();
    vertices = buffer;
}

// Difference with
void vertices_vec::setminus_with(vertices_vec const& v)
{
    // std::cout << fast_contain(7) << std::endl;
    // std::cout << v.fast_contain(7) << std::endl;

    // Put difference in the buffer
    buffer.clear();
    std::set_difference(vertices.begin(), vertices.end(), v.vertices.begin(),
        v.vertices.end(), std::inserter(buffer, buffer.begin()));
    // Copy the buffer in vertices
    // vertices.clear();
    vertices = buffer;

    // std::cout << fast_contain(7) << std::endl;
}

void vertices_vec::safe_difference_with(vertices_vec & v)
{
    if (!v.is_sorted())
        v.sort_by_number();
    if (!is_sorted())
        sort_by_number();
    difference_with(v);
}

/*
vertices_vec vertices_vec::union_with(vertices_vec const& v)
{
        vertices_vec union_vec(this);
        union_vec.union_with(v);
        return union_vec;
}

vertices_vec vertices_vec::safe_union_with(vertices_vec & v)
{
        vertices_vec union_vec(this);
        union_vec.safe_union_with(v);
        return union_vec;
}

vertices_vec vertices_vec::union_with(int i)
{
        vertices_vec union_vec(this);
        union_vec.union_with(i);
        return union_vec;
}
*/

void vertices_vec::add_interval(const int l, const int u) {
    vertices.clear();
    for (int i = l; i <= u; ++i) {
        vertices.push_back(i);
    }
}

std::ostream& vertices_vec::display(std::ostream& os) const {
    os << "{ ";
    for (std::vector<int>::const_iterator it = vertices.begin();
         it != vertices.end(); ++it) {
        os << *it << " ";
    }
    os << "}";
                return os;
}

void vertices_vec::remove(const int e) { vertices.erase(find(e)); }

bool vertices_vec::includes(vertices_vec const& v) const {
    return std::includes(begin(), end(), v.begin(), v.end());
}

bool vertices_vec::intersect(vertices_vec const& v) const {
    size_t i = 0;
    size_t j = 0;
    while (i < size() && j < v.size()) {
        if (vertices[i] == v.vertices[j])
            return true;
        if (vertices[i] > v.vertices[j])
            ++j;
        else
            ++i;
    }
    return false;
}

int vertices_vec::min() const { return *begin(); }

int vertices_vec::max() const { return *crbegin(); }

/********************************************************/
//                                      NEIGHBORS_VEC //
/********************************************************/

int neighbors_vec::degree() { return vertices.size(); }
} // namespace gc


std::ostream& gc::operator<<(std::ostream& os, const gc::vertices_vec& x)
{
    return x.display(os);
}


