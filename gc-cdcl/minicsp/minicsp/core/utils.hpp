/*************************************************************************
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

*************************************************************************/

#ifndef __MINICSP_UTILS_HPP__
#define __MINICSP_UTILS_HPP__

#include <stdint.h>

namespace minicsp {

class Solver;

double cpuTime(void);
uint64_t memUsed();
void printStats(Solver& solver, const char *comment = 0L);
void setup_signal_handlers(Solver *s);

}

#endif

