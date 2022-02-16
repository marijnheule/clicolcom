#ifndef __GC_STATISTICS_HPP
#define __GC_STATISTICS_HPP

#include <iostream>
#include "prop.hpp"
// #include <minicsp/core/solver.hpp>


using namespace minicsp;


namespace gc
{

void process_mem_usage(double& vm_usage, double& resident_set);

struct statistics {

    statistics(const int size)
    {

        changed = true;
        cons = NULL;

        // total_time = 0;
        total_conflicts = 0;
				total_iteration = 0;
        best_lb = 0;
        best_ub = size;

        update_lb = true;
        update_ub = true;

        ub_safe = true;

        total_bound_1 = 0;
        total_bound_2 = 0;

        total_clq_size = 0;
        num_bounds = 0;

        num_vertices = size;

        num_iterations = 0;
        core_size = 0;

        start_time = minicsp::cpuTime();
    }

    // outputs a nice description of all statistics
    void describe(std::ostream&);
    void display(std::ostream&);

    void binds(cons_base* c);
    void unbinds();

    cons_base* cons;

    bool changed;

    int best_lb;
    void notify_lb(const int l);
    int best_ub;
    void notify_ub(const int u);

    bool update_lb;
    bool update_ub;

    bool ub_safe;

    uint64_t total_clq_size;
    uint64_t num_bounds;
    int get_avg_nclq();
    void notify_nclique(const int sz);

    // double total_time;
    uint64_t total_conflicts;
		uint64_t total_iteration;

    // the actual statistics
    uint64_t num_neighborhood_contractions;
    // uint64_t num_vertex_removals;
    int num_vertices;
    void notify_removals(const int n);

    uint64_t total_bound_1;
    uint64_t total_bound_2;
    void notify_bound_delta(const int b1, const int b2);
    double get_bound_increase() const;

    uint64_t num_iterations;
    int core_size;
    void notify_iteration(const int cs);

    // in order to ignor reading time
    double start_time;
};

} // namespace gc

#endif
