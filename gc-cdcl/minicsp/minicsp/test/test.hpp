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

#ifndef __MINICSP_TEST
#define __MINICSP_TEST

#include "minicsp/core/solver.hpp"
#include <vector>
#include <map>
#include <string>
#include <stdexcept>

using namespace minicsp;

void assert_clause_exact(Solver &s, Clause *to_test,
                         vec<Lit> const& expected);
void assert_clause_contains(Solver &s, Clause *to_test,
                            vec<Lit> const& expected);
void assert_clause_exact(Solver &s, vec<Lit> &to_test,
                         vec<Lit> const& expected);
void assert_clause_contains(Solver &s, vec<Lit> &to_test,
                            vec<Lit> const& expected);

// the number of solutions of the model posted to s is exactly ns
void assert_num_solutions(Solver &s, int ns);

#define MUST_BE_UNSAT(x) do {                           \
    bool BOOST_PP_CAT(thrown, __LINE__) = false;        \
    try { x; } catch( unsat& ) {                        \
      BOOST_PP_CAT(thrown, __LINE__) = true; }          \
    assert(BOOST_PP_CAT(thrown, __LINE__));             \
  } while(0)

typedef void (*test)();

struct duplicate_test : public std::exception
{
  std::string tname;
  duplicate_test(std::string s) :tname(s) {}
  ~duplicate_test() throw() {}
  const char* what() const throw();
};

struct test_container {
  std::vector<std::pair<std::string, test> > tests;
  std::map< test, std::string> tset;
  void add(test t, std::string s);
  void run();
};

namespace { // diff singleton for each translation unit
  test_container& the_test_container() {
    static test_container t;
    return t;
  }

  struct register_actor {
    register_actor(test t, std::string s) {
      the_test_container().add(t, s);
    }
  };

}

#define REGISTER_TEST(x)                                 \
  register_actor BOOST_PP_CAT(register_actor_, __LINE__) \
    (x, #x)

#endif
