
// #include <unistd.h>
// #include <ios>
// #include <iostream>
#include <fstream>
// #include <string>
#include <iomanip>

#include "statistics.hpp"
#include <minicsp/core/utils.hpp>


namespace gc
{
	
	
	//////////////////////////////////////////////////////////////////////////////
	//
	// process_mem_usage(double &, double &) - takes two doubles by reference,
	// attempts to read the system-dependent data for a process' virtual memory
	// size and resident set size, and return the results in KB.
	//
	// On failure, returns 0.0, 0.0

	void process_mem_usage(double& vm_usage, double& resident_set)
	{
	   using std::ios_base;
	   using std::ifstream;
	   using std::string;

	   vm_usage     = 0.0;
	   resident_set = 0.0;

	   // 'file' stat seems to give the most reliable results
	   //
		 
		 
	   std::ifstream stat_stream("/proc/self/stat",std::ios_base::in);

		 if(stat_stream) {
	   // dummy vars for leading entries in stat that we don't care about
	   //
	   std::string pid, comm, state, ppid, pgrp, session, tty_nr;
	   std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	   std::string utime, stime, cutime, cstime, priority, nice;
	   std::string O, itrealvalue, starttime;

	   // the two fields we want
	   //
	   unsigned long vsize;
	   long rss;

	   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	               >> utime >> stime >> cutime >> cstime >> priority >> nice
	               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	   stat_stream.close();

	   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	   vm_usage     = vsize / 1024.0;
	   resident_set = rss * page_size_kb;
	 }
	}

	
		
void statistics::notify_bound_delta(const int b1, const int b2)
{
		total_bound_1 += b1;
		total_bound_2 += b2;
}

int statistics::get_avg_nclq()
{
    return (
        num_bounds ? (int)((double)total_clq_size / (double)num_bounds) : 0);
}

void statistics::notify_nclique(const int sz)
{
    total_clq_size += sz;
    ++num_bounds;
}

void statistics::notify_lb(const int l) 
{
    if (best_lb < l) {
        best_lb = l;
        changed = true;
    }
}

void statistics::notify_ub(const int u) 
{
    if (best_ub > u) {
        best_ub = u;
        changed = true;
    }
}

void statistics::notify_removals(const int n)
{
		if(num_vertices > n) {
				num_vertices = n;
				changed = true;
		}
}

double statistics::get_bound_increase() const {
		if(total_bound_1)
				return (double)total_bound_2/(double)total_bound_1;
		return 0;
}

void statistics::notify_iteration(const int cs)
{
    core_size = cs;
    ++num_iterations;
    if (num_iterations < 100 or num_iterations % 100 == 0) {
        changed = true;
        display(std::cout);
    }
}

void statistics::describe(std::ostream& os)
{
    os << "[statistics] lb ub clq iter core time conflicts moves memory\n";
}

void statistics::display(std::ostream& os)
{
    if (update_lb and cons and best_lb < cons->bestlb) {
        changed = true;
        best_lb = cons->bestlb;
    }
    if (ub_safe and update_ub and cons and best_ub > cons->ub) {
        changed = true;
        best_ub = cons->ub;
    }

    if (changed) {

        double vm_usage;
        double resident_set;

        process_mem_usage(vm_usage, resident_set);

        os.setf(std::ios_base::fixed, std::ios_base::floatfield);
        os << "[data] lb = " << std::setw(4) << std::left << best_lb
           << "| ub = " << std::setw(4) << std::left << best_ub
           << "| clq = " << std::setw(4) << std::left << get_avg_nclq()
           << "| iter = " << std::setw(4) << std::left << num_iterations
           << "| core = " << std::setw(4) << std::left << core_size
           << "| time = " << std::setw(9) << std::left << std::setprecision(4)
           << (minicsp::cpuTime() - start_time)
           << "| conflicts = " << std::setw(8) << std::left
           << (cons ? total_conflicts + cons->s.conflicts : total_conflicts)
           << "| moves = " << std::setw(10) << total_iteration
           << "| memory = " << std::setw(8) << std::left << (long)resident_set
           << std::endl;
    }

    changed = false;
}

void statistics::binds( gc::cons_base* c ) {
		cons = c;
}

void statistics::unbinds() {
		if(cons)
				total_conflicts += cons->s.conflicts;
		cons = NULL;
}

} // namespace gc
