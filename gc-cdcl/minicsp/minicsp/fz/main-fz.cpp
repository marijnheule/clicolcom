/*****************************************************************************
  minicsp

Copyright 2010 George Katsirelos

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

*****************************************************************************/

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <iomanip>

#include "flatzinc.hpp"
#include "minicsp/core/solver.hpp"
#include "minicsp/core/cmdline.hpp"
#include "minicsp/core/utils.hpp"

using namespace std;
using namespace minicsp;

int main(int argc, char *argv[])
{
  list<string> args(argv+1, argv+argc);

  if(args.empty()) {
    cout << "usage: minicsp-fz [options] input.fzn";
    return 1;
  }

  Solver s;

  cmdline::parse_solver_options(s, args);
  bool stat = cmdline::has_option(args, "--stat");
  bool maint = cmdline::has_option(args, "--maint");
  setup_signal_handlers(&s);

  double cpu_time = cpuTime();

  FlatZinc::Printer p;
  FlatZinc::FlatZincModel *fm = 0L;
  if( maint ) // do not catch exceptions
    fm = parse(args.back(), s, p);
  else try {
      fm = parse(args.back(), s, p);
    } catch (unsat& e) {
      cout << setw(5) << setfill('=') << '='
           << "UNSATISFIABLE" << setw(5) << '=' << "\n";
    }
  if( !fm ) return 0;

  double parse_time = cpuTime() - cpu_time;

  fm->findall = cmdline::has_option(args, "--all");

  fm->run(cout , p);
  delete fm;

  if( stat ) {
    printStats(s, "%% ");
    reportf("%sParse time            : %g s\n", "%% ", parse_time);
  }

  return 0;
}
