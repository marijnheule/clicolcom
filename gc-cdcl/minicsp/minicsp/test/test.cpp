/*************************************************************************
minicsp

Copyright 2010--2011 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/

#include "test.hpp"
#include "minicsp/core/solver.hpp"

using namespace minicsp;

template<typename T>
void assert_clause_exact0(Solver &s,
                          T& to_test, vec<Lit> const& expected)
{
  if( expected.size() != to_test.size() )
    std::cout << "Expected " << print(s, &expected)
              << "\ngot " << print(s, &to_test) << "\n";
  assert(expected.size() == to_test.size());
  assert_clause_contains0(s, to_test, expected);
}

template<typename T>
void assert_clause_contains0(Solver &s,
                             T& to_test, vec<Lit> const& expected)
{
  for(int i = 0; i != expected.size(); ++i) {
    bool f = false;
    for(int j = 0; j != to_test.size(); ++j)
      if( expected[i] == to_test[j] ) {
        f = true;
        break;
      }
    if( !f )
      std::cout << "Expected " << print(s, &expected)
                << "\ngot " << print(s, &to_test) << "\n";
    assert(f);
  }
}


void assert_clause_exact(Solver &s, Clause *to_test,
                         vec<Lit> const& expected)
{
  assert_clause_exact0(s, *to_test, expected);
}

void assert_clause_contains(Solver &s, Clause *to_test,
                            vec<Lit> const& expected)
{
  assert_clause_contains0(s, *to_test, expected);
}

void assert_clause_exact(Solver &s, vec<Lit> &to_test,
                         vec<Lit> const& expected)
{
  assert_clause_exact0(s, to_test, expected);
}

void assert_clause_contains(Solver &s, vec<Lit> &to_test,
                            vec<Lit> const& expected)
{
  assert_clause_contains0(s, to_test, expected);
}

void assert_num_solutions(Solver &s, int ns)
{
  int numsol = 0;
  bool next;
  do {
    next = s.solve();
    if( next ) {
      ++numsol;
      try {
        s.excludeLast();
      } catch(unsat&) {
        next = false;
      }
    }
  } while(next);
  assert(numsol == ns);
}

const char *duplicate_test::what() const throw()
{
  return tname.c_str();
}

void test_container::add(test t, std::string s)
{
  if( tset.find(t) != tset.end() )
    throw duplicate_test( tset[t]);
  tests.push_back(make_pair(s, t));
  tset[t] = s;
}

void test_container::run()
{
  using namespace std;
  for(size_t i = 0; i != tests.size(); ++i) {
    cout << tests[i].first << "..." << flush;
    (*tests[i].second)();
    cout << "OK" << endl;
  }
}
