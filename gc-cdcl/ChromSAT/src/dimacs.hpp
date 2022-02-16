#ifndef __CG_DIMACS_HH
#define __CG_DIMACS_HH

#include <fstream>
#include <sstream>

namespace dimacs
{

template <typename setsize, typename add_edge, typename add_node>
void read_graph(
    const char* fn, setsize ss, add_edge ae, add_node an, bool strict = false)
{
    using std::cerr;
    try {
        std::ifstream ifs(fn);
        if (!ifs)
            throw std::runtime_error("Could not open file for reading");
        bool gotheader{false};
        std::string lt; // line type
        int ln{1};
        int nv{0}, ne{0}, reade{0};
        int64_t w;
        int checkcount = 0;
        for (std::string line; getline(ifs, line); ++ln) {
            if (line[0] != 'p' && line[0] != 'e')
                continue;
            std::istringstream iss(line);
            if (!gotheader) {
                std::string edge;
                iss >> lt >> edge >> nv >> ne;
                if (!iss || lt != "p" || (edge != "edge" && edge != "edges" && edge != "col")) {
                    cerr << "ERROR: could not parse header at line " << ln
                         << "\n";
                    exit(1);
                }
                ss(nv, ne);
                gotheader = true;
            } else {
                // edge
                int v1, v2;
                iss >> lt >> v1 >> w;
                if (lt == "n") {
                    ++checkcount;

                    if (checkcount != v1) {
                        cerr << "ERROR: wrong node numbering at line " << ln
                             << "\n";
                        exit(1);
                    }

                    an(v1, w);

                } else {

                    if (lt != "e") {
                        cerr << "ERROR: could not parse edge at line " << ln
                             << "\n";
                        exit(1);
                    }
                    v2 = w;

                    if (strict && (v1 > nv || v2 > nv)) {
                        cerr << "ERROR: " << nv
                             << " vertices declared, but edge " << v1 << " -- "
                             << v2 << " requires more at line " << ln << "\n";
                        exit(1);
                    }
                    ++reade;
                    if (strict && reade > ne) {
                        cerr << "ERROR: " << ne
                             << " edges declared, but at least " << reade
                             << " edges present at line " << ln << "\n";
                        exit(1);
                    }
										
										//
										// if(v1<1 or v2<1 or v1>nv or v2>nv)
										// {
										// 	std::cout << "edge " << v1 << "." << v2 << " at line " << ln << std::endl;
										// 	exit(1);
										// }
										
                    ae(v1, v2);
                }
            }
        }
        if (strict && reade < ne) {
            cerr << "ERROR: " << ne << " edges declared, but only " << reade
                 << " edges read\n";
            exit(1);
        }
    } catch (std::exception& e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
} // namespace dimacs

#endif
