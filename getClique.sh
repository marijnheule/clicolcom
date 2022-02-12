
GRAPH=$1

./cliquer/cl $GRAPH | grep "=" | tr " " "\n" | grep -v "=" | awk '{if (NF == 1) print $0}'
