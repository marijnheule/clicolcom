
#include <iostream>

#include <vector>

#ifndef __INTSTACK_HPP
#define __INTSTACK_HPP

/**********************************************
 * intstack
 **********************************************/
/// Sparse set representation

class intstack
{
private:
    /*!@name Parameters*/
    //@{
    /// list of values
    std::vector<int> list_;
    size_t size_;

    /// values' indices
    std::vector<size_t> index_;
    //@}

public:
    /*!@name Constructors*/
    //@{
    explicit intstack(const size_t n = 0);
    explicit intstack(std::vector<size_t>& common);

    void reserve(const size_t n);

    /*!@name Accessors*/
    //@{
    bool safe_contain(const int elt) const;
    bool contain(const int elt) const;

    size_t size() const;
    bool empty() const;

    int next(const int elt) const;

    int prev(const int elt) const;

    int operator[](const size_t idx) const;

    int& operator[](const size_t idx);
    //@}

    /*!@name List Manipulation*/
    //@{
    std::vector<int>::iterator begin();
    std::vector<int>::reverse_iterator rbegin();

    std::vector<int>::iterator end();
    std::vector<int>::reverse_iterator rend();

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_reverse_iterator rbegin() const;

    std::vector<int>::const_iterator end() const;
    std::vector<int>::const_reverse_iterator rend() const;

    std::vector<int>::iterator begin_not_in();
    std::vector<int>::reverse_iterator rbegin_not_in();

    std::vector<int>::iterator end_not_in();
    std::vector<int>::reverse_iterator rend_not_in();

    std::vector<int>::const_iterator begin_not_in() const;
    std::vector<int>::const_reverse_iterator rbegin_not_in() const;

    std::vector<int>::const_iterator end_not_in() const;
    std::vector<int>::const_reverse_iterator rend_not_in() const;

    void fill();

    void clear();

    void resize(const size_t n);

    void move_up(const int elt, const int idx);

    void pop_back();

    void pop_head();

    int head() const;

    int back() const;

    void push(const int elt);
    void add(const int elt);
    void safe_add(const int elt);

    void pull(const int elt);
    void remove(const int elt);
    void safe_remove(const int elt);

    int index(const int elt) const;

    void save(size_t& stamp);
    void restore(const size_t stamp);
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const intstack& x);

#endif // __INTSTACK_HPP
