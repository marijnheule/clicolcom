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

#include "XCSP3CoreParser.h"
#include "XCSP3MiniCSPCallbacks.hpp"

using namespace XCSP3Core;
using namespace std;
using namespace minicsp;

// helper for variant, from cppreference
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

int main(int argc, char *argv[]) {
    list<string> args(argv + 1, argv + argc);

    if(args.empty()) {
        cout << "usage: minicsp-xcsp3 [options] input.xml";
        return 1;
    }

    Solver s;

    cmdline::parse_solver_options(s, args);
    setup_signal_handlers(&s);

    double cpu_time = cpuTime();
    XCSP3MiniCSPCallbacks cb(s); // my interface between the parser and the solver

    pair<bool, string> has_removeClasses = cmdline::has_argoption<string>(args, "--rmclass");
    if(has_removeClasses.first) {
        std::vector<std::string> classes = split(std::string(has_removeClasses.second), ',');
        for(std::string c : classes)
            cb.addClassToDiscard(c);
    }
    try {
        XCSP3CoreParser parser(&cb);
        parser.parse(args.back().c_str()); // fileName is a string
    }
    catch(exception &e) {
        cout.flush();
        cerr << "\n\tUnexpected exception :\n";
        cerr << "\t" << e.what() << endl;
        exit(1);
    }


    double parse_time = cpuTime() - cpu_time;
    cout << parse_time << " s to parse instance" << endl;

    bool findall = cmdline::has_option(args, "--all");

    bool optimization = !holds_alternative<ConstantObjective>(cb.objective);

    bool sat = false, next;
    bool optimum = false;
    int ns = 0;
    do {
      next = s.solve();
      if (!next)
          optimum = true;
      sat = sat || next;
      if (next) {
        ++ns;
        cout << "c solution " << ns << ": ";
        cb.print_solution();
        try {
          visit(overloaded{[&](ConstantObjective) {
                             if (findall) {
                               s.excludeLast();
                               next = findall;
                             }
                           },
                           [&](MinimizeObjective obj) {
                             int curval = s.cspModelValue(obj.var);
                             cout << "o " << curval << "\n";
                             obj.var.setmax(s, curval - 1, NO_REASON);
                             next = true;
                           },
                           [&](MaximizeObjective obj) {
                             int curval = s.cspModelValue(obj.var);
                             cout << "o " << curval << "\n";
                             obj.var.setmin(s, curval + 1, NO_REASON);
                             next = true;
                           }},
                cb.objective);
        } catch (unsat &e) {
          next = false;
          optimum = true;
        }
      }
    } while (next);

    if (!optimization) {
      if (sat)
        cout << "s SATISFIABLE\n";
      else
        cout << "s UNSATISFIABLE\n";
    } else {
      if (sat && optimum)
        cout << "s OPTIMUM FOUND\n";
      else if (sat && !optimum)
        cout << "s SATISFIABLE\n";
      else if (!sat)
        cout << "s UNSATISFIABLE\n";
    }

    printStats(s, "c ");
    reportf("%sParse time            : %g s\n", "c ", parse_time);

    return 0;
}
