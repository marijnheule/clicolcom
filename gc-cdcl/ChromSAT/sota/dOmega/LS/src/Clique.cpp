/**
 * @file Clique.cpp
 *
 * @brief Finds The maximum clique of a graph.
 *
 * @details This code implements an efficient algorithm for the maximum clique
 * problem that runs in time polynomial in the graph's size, but exponential
 * in the difference $d-\omega$ (where $d$ denotes the graph's degeneracy).
 * When this difference is a constant, the running time is ${O}(dm)={O}(m^{1.5})$.
 * Since 70% of common test instances have $d-\omega \le 3$, our implementation
 * performs rather well and is competitive with previous approaches.
 * Key subroutines in our approach include existing kernelization and
 * fixed-parameter tractable algorithms for vertex cover.
 *
 * Please cite the following paper if used:
 * Jose L. Walteros and Austin L. Buchanan. Why is maximum clique easy in practice?
 *
 * @author Jose L. Walteros
 *
 * @version   1.0
 * @date      June 2018
 *
 * Copyright (C) 2018 Jose L. Walteros. All rights reserved.
 *
 *---+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

#include <math.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <algorithm>
#include "Clique.h"
#include "Graph.h"
#include "NemhauserTrotter.h"
#include "Buss.h"
#include "VertexCover.h"

namespace minicsp
{
double cpuTime();
}
using minicsp::cpuTime;

using namespace dOmega;

void processSubgraphs(Graph& graph, std::vector<int>& sortedList,
    std::vector<subgraph>& subgraphs, std::atomic<bool>& cliqueFlag,
    std::atomic<bool>& interrupted, int threadNumber, int numberOfThreads,
    int clq, double timeout = -1)
{
    VertexCover VC;
    int i = threadNumber;

    while (!cliqueFlag && i < graph.n) {
        int v = sortedList[i];

        int k = graph.rightDegree[v] + 1 - clq;

        if (k >= 0) {
            /**
             * If the subgraph of vertex v has not been created, it does it.
             */
            if (subgraphs[v].created == false) {
                if (timeout > 0 && cpuTime() > timeout) {
                    interrupted = true;
                    return;
                }
                graph.generateCompGraphRightNeighbors(v, subgraphs);
            }

            /**
             * Generates the Buss kernel.
             */
            Buss BusKernel(&subgraphs[v], k);
            subgraph kernel;
            int highDegVertices = 0;
            int success = BusKernel.getKernel(kernel, highDegVertices);

            if (success == -1) {
                i += numberOfThreads;
                continue;
            }

            if (success == 1) {
                cliqueFlag = true;
                break;
            }
            k = k - highDegVertices;

            /**
             * Generates the NT kernel.
             */

            subgraph kernel2;
            int numRemoved = 0;
            int numInVC = 0;
            NemhauserTrotter NT(&kernel, k);
            success = NT.getKernel(kernel2, numRemoved, numInVC);

            if (success == -1) {
                i += numberOfThreads;
                continue;
            }

            if (success == 1) {
                cliqueFlag = true;
                break;
            }

            k = k - numInVC;

            if (timeout > 0 && cpuTime() > timeout) {
                interrupted = true;
                return;
            }

            /**
             * Solves the resulting k vertex cover problem.
             */
            if (VC.kVertexCover(kernel2.n, k, kernel2.vertices,
                    kernel2.adjLists, interrupted, timeout)) {
                cliqueFlag = true;
                break;
            }
        } else {
            break;
        }
        i += numberOfThreads;
    }
}

Clique::Clique(
    Graph& graph,
    const int numThreads) : graph(graph) {
    this->numThreads = numThreads;
}

int Clique::findMaxClique(int lb, int ub, double timeout)
{
    subgraphs = std::vector<subgraph>(graph.n);
    std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();
    graph.degeneracyOrdering(subgraphs);
    cliqueUB = std::min(ub, graph.cliqueUB);
    cliqueLB = std::max(lb, graph.cliqueLB);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    degeneracyTime = std::chrono::duration_cast<std::chrono::duration<double> >(end_time - begin_time);

    /**
     * If the upper and lower bounds are different, the vertices are sorted
     * based on their right degree
     */

    if (cliqueLB < cliqueUB) {
        sortedList = std::vector<int>(graph.n);
        std::vector<int>buckets(graph.d + 1, 0);

        for (int i = 0; i < graph.n; i++) {
            buckets[graph.rightDegree[i]]++;
        }

        int temp;
        int count = 0;

        for (int k = graph.d; k >= 0; k--) {
            temp       = buckets[k];
            buckets[k] = count;
            count     += temp;
        }

        for (int i = 0; i < graph.n; i++) {
            sortedList[buckets[graph.rightDegree[i]]] = i;
            buckets[graph.rightDegree[i]]++;
        }

        int clq = cliqueUB;
        std::vector<std::thread>threads(numThreads);

        while (cliqueLB < cliqueUB) {
            std::cout << "[info] maxclique in [" << cliqueLB << "," << cliqueUB
                      << "]\n";
            cliqueFlag = false;
            interrupted = false;

            for (int i = 0; i < numThreads; i++) {
                std::thread th(&processSubgraphs, std::ref(graph),
                    std::ref(sortedList), std::ref(subgraphs),
                    std::ref(cliqueFlag), std::ref(interrupted), i, numThreads,
                    clq, timeout);
                threads[i] = std::move(th);
            }


            for (std::thread& th : threads) {
                th.join();
            }

            if (interrupted) {
                break;
            }

            if (cliqueFlag) {
                cliqueLB = clq;
            } else {
                cliqueUB = clq - 1;
            }
            clq = cliqueUB;
        }
    }
    // end_time  = std::chrono::high_resolution_clock::now();
    // runningTime = std::chrono::duration_cast<std::chrono::duration<double> >(end_time - begin_time);

    // std::clog << "Number of threads used: " << numThreads << "\n";
    // std::clog << "Degeneracy: " << graph.d << "\n";
    // std::clog << "Lower bound from degeneracy: " << graph.cliqueLB << "\n";
    // std::clog << "Maximum clique size: " << cliqueUB << "\n";
    // std::clog << "Total running time: " << runningTime.count() << " \n";
    // std::clog << "-------------------------------------------------------------\n";

    return 0;
}
