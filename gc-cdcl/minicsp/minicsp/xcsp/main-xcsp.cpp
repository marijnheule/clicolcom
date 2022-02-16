/******************************************************************************
minicsp

Copyright 2011 George Katsirelos

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

Parts of this file were distributed under the following license as
part of the xcsp distribution

 Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
******************************************************************************/

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "minicsp/core/solver.hpp"
#include "minicsp/core/cmdline.hpp"
#include "minicsp/core/utils.hpp"

#include "XMLParser_libxml2.hh"
#include "minicsp_callback.hpp"

using namespace CSPXMLParser;

using namespace std;

int main(int argc, char *argv[])
{
  list<string> args(argv+1, argv+argc);

  if(args.empty()) {
    cout << "usage: minicsp-xcsp [options] input.xml";
    return 1;
  }

  Solver s;

  cmdline::parse_solver_options(s, args);
  bool stat = cmdline::has_option(args, "--stat");
  bool maint = cmdline::has_option(args, "--maint");
  setup_signal_handlers(&s);

  double cpu_time = cpuTime();

  minicsp_callback cb(s);

  XMLParser_libxml2<> parser(cb);
  parser.setPreferredExpressionRepresentation(TREE);
  if( maint ) // do not catch exceptions
    parser.parse(args.back().c_str());
  else try {
      parser.parse(args.back().c_str());
    } catch (unsat& e) {
      cout << setw(5) << setfill('=') << '='
           << "UNSATISFIABLE" << setw(5) << '=' << "\n";
      return 0;
    } catch (unsupported& e) {
      cout << "\nERROR: unsupported instance\n";
      return 1;
    }
  double parse_time = cpuTime() - cpu_time;
  cout << parse_time << " s to parse instance\n";

  bool findall = cmdline::has_option(args, "--all");

  bool sat = false, next;
  int ns = 0;
  do {
    next = s.solve();
    sat = sat || next;
    if( next ) {
      ++ns;
      cout << "solution " << ns << ": ";
      cb.print_solution();
      if ( findall ) {
        try {
          s.excludeLast();
        } catch (unsat& e) {
          next = false;
        }
      }
    }
    next = next && findall;
  } while(next);

  if( sat )
    cout << "sat\n";
  else
    cout << "unsat\n";

  if( stat ) {
    printStats(s, "%% ");
    reportf("%sParse time            : %g s\n", "%% ", parse_time);
  }

  return 0;
}
