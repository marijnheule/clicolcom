#ifndef __CG_EDGEFORMAT_HH
#define __CG_EDGEFORMAT_HH

#include <fstream>
#include <sstream>

namespace edgeformat
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
        int nv{0}, ne{0}, u, v;
				ifs >> nv >> ne;
				ss(nv, ne);
				
				for(int i=0; i<ne; ++i) {
					ifs >> u >> v;
					ae(u, v);
				}
    } catch (std::exception& e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
} // namespace edgeformat

#endif
