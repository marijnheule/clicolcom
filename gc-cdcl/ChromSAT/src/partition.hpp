
#include <iostream>

#include <vector>

#ifndef __PARTITION_HPP
#define __PARTITION_HPP

/**********************************************
 * partition
 **********************************************/
/// Sparse set representation

#define NOTIN 0xfffffff

class partition
{
private:
  /// values' indices
  std::vector<size_t> index_;
  //@}
	
public:
    /*!@name Parameters*/
    //@{
    /// list of values
    std::vector<std::vector<int>> bag;

    /*!@name Constructors*/
    //@{
    explicit partition();

    void clear();
    void resize(const size_t n, const size_t m);

    int size();

    const std::vector<int>& operator[](const int i) const;

    /*!@name List Manipulation*/
    //@{
    void move(const int elt, const int from, const int to);
    void add_elt(const int elt, const int to);
    void swap(const int a, const int b);
    void remove(const int a);
    bool contain(const int elt, const int c);
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const partition& x);

#endif // __PARTITION_HPP
