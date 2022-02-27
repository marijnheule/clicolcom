GRAPH=$1
timeout 1 ./inc_max_clique $GRAPH | grep "Clique Output:" | tail -n 1 | sed 's|Clique Output: ||' | tr " " "\n"
