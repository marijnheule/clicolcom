/******************************************************************************************[Main.C]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <ctime>
#include <cstring>
#include <stdint.h>
#include <errno.h>

#include <string>
#include <iostream>

#include <signal.h>
#include <zlib.h>

#include "Solver.h"

/*************************************************************************************/
#ifdef _MSC_VER
#include <ctime>

static inline double cpuTime(void) {
    return (double)clock() / CLOCKS_PER_SEC; }
#else

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

static inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }
#endif


#if defined(__linux__)
static inline int memReadStat(int field)
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
static inline uint64_t memUsed() { return (uint64_t)memReadStat(0) * (uint64_t)getpagesize(); }


#elif defined(__FreeBSD__)
static inline uint64_t memUsed(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss*1024; }


#else
static inline uint64_t memUsed() { return 0; }
#endif

#if defined(__linux__)
#include <fpu_control.h>
#endif

//=================================================================================================
// DIMACS Parser:

#define CHUNK_LIMIT 1048576

char tmp_file[512];

class StreamBuffer {
    gzFile  in;
    char    buf[CHUNK_LIMIT];
    int     pos;
    int     size;

    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = gzread(in, buf, sizeof(buf)); } }

public:
    StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
        assureLookahead(); }

    int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
    void operator ++ () { pos++; assureLookahead(); }
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<class B>
static void skipWhitespace(B& in) {
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; }

template<class B>
static void skipLine(B& in) {
    for (;;){
        if (*in == EOF || *in == '\0') return;
        if (*in == '\n') { ++in; return; }
        ++in; } }

template<class B>
static int parseInt(B& in) {
    int     val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }

template<class B>
static void readClause(B& in, Solver& S, vec<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? Lit(var) : ~Lit(var) );
    }
}

template<class B>
static bool match(B& in, char* str) {
    for (; *str != 0; ++str, ++in)
        if (*str != *in)
            return false;
    return true;
}


template<class B>
static void parse_DIMACS_main(B& in, Solver& S) {
    vec<Lit> lits;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF)
            break;
        else if (*in == 'p'){
            if (match(in, "p cnf")){
                int vars    = parseInt(in);
                int clauses = parseInt(in);
                // reportf("|  Number of variables:  %-12d                                         |\n", vars);
                // reportf("|  Number of clauses:    %-12d                                         |\n", clauses);
            }else{
                reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }

        } else if (*in == 'e') {
          if (match(in, "e info")) {
            int v = parseInt(in);
            int c = parseInt(in);
            S.initVertexVars(v,c);
            // reportf("|  Number of nodes:      %-12d                                         |\n", v);
            // reportf("|  Number of colors:     %-12d                                         |\n", c);
          }
        } else if (*in == 'c' || *in == 'p') {
            skipLine(in);
        } else {
            readClause(in, S, lits),
            S.addClause(lits);
        }
    }
}

// Inserts problem into solver.
//
static void parse_DIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S); }


//=================================================================================================


void printStats(Solver& solver)
{
    double   cpu_time = cpuTime();
    uint64_t mem_used = memUsed();
    reportf("restarts              : %lld\n", solver.starts);
    reportf("conflicts             : %-12lld   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    reportf("decisions             : %-12lld   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    reportf("propagations          : %-12lld   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    reportf("conflict literals     : %-12lld   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) reportf("Memory used           : %.2f MB\n", mem_used / 1048576.0);
    reportf("CPU time              : %g s\n", cpu_time);
}

Solver* solver;
static void SIGINT_handler(int signum) {
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
		// system((char *)tmp_file);
    printStats(*solver);
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    exit(1); }


//=================================================================================================
// Main:

void printUsage(char** argv)
{
    reportf("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n\n", argv[0]);
    reportf("OPTIONS:\n\n");
    reportf("  -polarity-mode = {true,false,rnd}\n");
    reportf("  -decay         = <num> [ 0 - 1 ]\n");
    reportf("  -rnd-freq      = <num> [ 0 - 1 ]\n");
    reportf("  -verbosity     = {0,1,2}\n");
    reportf("\n");
}


const char* hasPrefix(const char* str, const char* prefix)
{
    int len = strlen(prefix);
    if (strncmp(str, prefix, len) == 0)
        return str + len;
    else
        return NULL;
}

bool solve( gzFile& in, int argc, char** argv, double& search_time, double& parse_time, long int& n_conflicts ) {
	
	  Solver      S;
	  S.verbosity = -1;


	  int         i, j;
	  const char* value;
	  for (i = j = 0; i < argc; i++){
	      if ((value = hasPrefix(argv[i], "-polarity-mode="))){
	          if (strcmp(value, "true") == 0)
	              S.polarity_mode = Solver::polarity_true;
	          else if (strcmp(value, "false") == 0)
	              S.polarity_mode = Solver::polarity_false;
	          else if (strcmp(value, "rnd") == 0)
	              S.polarity_mode = Solver::polarity_rnd;
	          else{
	              reportf("ERROR! unknown polarity-mode %s\n", value);
	              exit(0); }

	      }else if ((value = hasPrefix(argv[i], "-rnd-freq="))){
	          double rnd;
	          if (sscanf(value, "%lf", &rnd) <= 0 || rnd < 0 || rnd > 1){
	              reportf("ERROR! illegal rnd-freq constant %s\n", value);
	              exit(0); }
	          S.random_var_freq = rnd;

	      }else if ((value = hasPrefix(argv[i], "-decay="))){
	          double decay;
	          if (sscanf(value, "%lf", &decay) <= 0 || decay <= 0 || decay > 1){
	              reportf("ERROR! illegal decay constant %s\n", value);
	              exit(0); }
	          S.var_decay = 1 / decay;

	      }else if ((value = hasPrefix(argv[i], "-verbosity="))){
	          int verbosity = (int)strtol(value, NULL, 10);
	          if (verbosity == 0 && errno == EINVAL){
	              reportf("ERROR! illegal verbosity level %s\n", value);
	              exit(0); }
	          S.verbosity = verbosity;

	      }else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0){
	          printUsage(argv);
	          exit(0);

	      }else if (strncmp(argv[i], "-", 1) == 0){
	          reportf("ERROR! unknown flag %s\n", argv[i]);
	          exit(0);

	      }else
	          argv[j++] = argv[i];
	  }
	  argc = j;
		
	
		double cpu_time = cpuTime();
	
	  solver = &S;
	  signal(SIGINT,SIGINT_handler);
	  signal(SIGHUP,SIGINT_handler);


	  parse_DIMACS(in, S);
	  gzclose(in);

		double time_now = cpuTime();
	  parse_time += time_now - cpu_time;


	  if (!S.isInitialized()) {
	    reportf("ERROR no vertex parameters were provided!\n");
			// return false;
	    // exit(0);
	  }
	  if (!S.simplify()){
	      reportf("Solved by unit propagation\n");
	      printf("UNSATISFIABLE\n");
				return false;
	      // exit(20);
	  }

	  bool ret = S.solve();
		
		search_time += (cpuTime() - time_now);
		
		n_conflicts += S.conflicts;
		
	  // printStats(S);
	  // reportf("\n");
	  printf(ret ? "SATISFIABLE\n" : "UNSATISFIABLE\n");

	
		return ret;
}


int main(int argc, char** argv)
{
		double search_time = 0;		
		double parse_time = 0;
		double encode_time = 0;
		
		long int n_conflicts = 0;
		
    if (argc != 5)
        reportf("Usage: %s <col file> <lb> <ub> <algo in {0 (TOP_DOWN), 1 (BOTTOM_UP), 2 (BINARY)}>\n");
			
		printf(">>statistics: lb ub time parsetime encodetime conflicts\n");
		
		time_t rawtime;
		struct tm * timeinfo;
		char temp[512];
		time (&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(temp,sizeof(temp),"%d-%m-%Y_%I:%M:%S",timeinfo);
		std::string timeofday(temp);
		std::string filename(argv[1]);
		size_t dot_pos = filename.find_last_of('.');
		std::string basename = filename.substr(0,dot_pos) + timeofday;
		std::string cnfname = basename + std::string(".cnf");	
		sprintf(tmp_file, "rm %s", cnfname.c_str());


		int lb = atoi(argv[2]);
		int ub = atoi(argv[3]);
		
		enum algo_type { TOP_DOWN, BOTTOM_UP, BINARY };
		
		
		int algo = atoi(argv[4]);
		algo_type policy{(algo == 0 ? TOP_DOWN : (algo == 1 ? BOTTOM_UP : BINARY))};
		
		double x_a = 0.0;
		double x_b = 0.0;
		double x_c = 0.0;
		long int x_d = 0;
		
		printf( ">>data: lb = %8d, ub = %8d, time = %10f, parsetime = %10f, encodetime = %10f, conflicts = %8ld\n", 
						lb, ub, x_a, x_b, x_c, x_d);
		
		
		bool solved = (lb == ub);
		
		if( ub < 32 ) {
				int n_colors = ( policy == TOP_DOWN ? ub-1 : ( policy == BOTTOM_UP ? lb : (lb+ub)/2 ) );
				while( lb < ub ) {
		
						printf("\nc try with #colors = %d (%s)\n", n_colors, basename.c_str());
		
						double cpu_time = cpuTime();
						sprintf(temp, "./sota/converter/converter %s %d 1 1 %s", filename.c_str(), n_colors, basename.c_str());
				
						// printf( "%s", temp );
				
						system((char *)temp);
						encode_time += cpuTime() - cpu_time;

				    gzFile in = gzopen(cnfname.c_str(), "rb");
				    if (in == NULL)
				        reportf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
		
						bool success = solve( in, argc, argv, search_time, parse_time, n_conflicts );
		
		
						if( success ) {
								ub = n_colors--;
						} else {
								lb = ++n_colors;
						}
						
						printf( ">>data: lb = %8d, ub = %8d, time = %10f, parsetime = %10f, encodetime = %10f, conflicts = %8ld\n", 
										lb, ub, search_time, parse_time, encode_time, n_conflicts);
				
						if(policy == BINARY) n_colors = (lb+ub)/2;
		
				}
		}
		
		
		printf( "\nbest bounds %d %d\n", lb, ub);

		// sprintf(temp, "rm %s", cnfname.c_str());
		
		// if( !solved )
		// 		system((char *)tmp_file);

}
