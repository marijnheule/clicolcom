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

#include "utils.hpp"
#include "solver.hpp"

#include <signal.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#ifdef _MSC_VER
#include <ctime>
#endif

namespace minicsp {

/*************************************************************************************/
#ifdef _MSC_VER
double cpuTime(void) {
    return (double)clock() / CLOCKS_PER_SEC; }
#else
double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }
#endif


#if defined(__linux__)
int memReadStat(int field)
{
    char    name[256];
    pid_t pid = getpid();
    sprintf(name, "/proc/%d/statm", pid);
    FILE*   in = fopen(name, "rb");
    if (in == NULL) return 0;
    int     value;
    for (; field >= 0; field--)
        fscanf(in, "%d", &value);
    fclose(in);
    return value;
}
uint64_t memUsed() { return (uint64_t)memReadStat(0) * (uint64_t)getpagesize(); }


#elif defined(__FreeBSD__)
uint64_t memUsed(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss*1024; }


#else
uint64_t memUsed() { return 0; }
#endif

void printStats(Solver& solver, const char *comment)
{
    if( !comment )
      comment = "";

    double   cpu_time = cpuTime();
    uint64_t mem_used = memUsed();
    reportf("%sBoolean variables     : %d\n", comment, solver.nVars());
    reportf("%sCSP variables         : %d\n", comment, solver.nCSPVars());
    reportf("%sClauses               : %d\n", comment, solver.nClauses());
    reportf("%sConstraints           : %d\n", comment, solver.nConstraints());
    reportf("%srestarts              : %" PRIu64 "\n", comment,
            solver.starts);
    reportf("%sconflicts             : %-12" PRIu64 "   (%.0f /sec)\n", comment,
            solver.conflicts   , solver.conflicts   /cpu_time);
    reportf("%sdecisions             : %-12" PRIu64 "   (%4.2f %% random) "
            "(%.0f /sec)\n",
            comment,
            solver.decisions,
            (float)solver.rnd_decisions*100 / (float)solver.decisions,
            solver.decisions   /cpu_time);
    reportf("%spropagations          : %-12" PRIu64 "   (%.0f /sec)\n",
            comment,
            solver.propagations, solver.propagations/cpu_time);
    reportf("%sconflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n",
            comment,
            solver.tot_literals,
            (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0)
      reportf("%sMemory used           : %.2f MB\n",
              comment, mem_used / 1048576.0);
    reportf("%sCPU time              : %g s\n", comment, cpu_time);
}

static Solver* global_solver;
static double previous_cputime;
#ifndef _MSC_VER
static void SIGINT_handler(int signum) {
    double curtime = cpuTime();
    if( curtime - previous_cputime < 2 ) {
      reportf("\n"); reportf("*** EXITING ***\n");
      exit(0);
    }
    previous_cputime = curtime;

    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    printStats(*global_solver);
}

static void SIGUSR_handler(int signum) {
    reportf("\n"); reportf("*** SIGUSR1 received ***\n");
    printStats(*global_solver);
}
#endif

void setup_signal_handlers(Solver *s)
{
  global_solver = s;
  previous_cputime = -5;
#ifndef _MSC_VER
  signal(SIGINT,SIGINT_handler);
  signal(SIGHUP,SIGINT_handler);
  signal(SIGXCPU,SIGINT_handler);
  signal(SIGUSR1,SIGUSR_handler);
#endif
}

} //namespace minicsp
