GRAPH=$1
SEED=$2
./inc_max_clique $GRAPH | grep Output | tail -n 1 | awk '{$1=""; $2=""; print $0}' |sed 's|^  ||' | tr " " "\n" > tmp-$$.clq
SIZE=`wc tmp-$$.clq | awk '{print $1}'`
./color $GRAPH $SIZE tmp-$$.clq | ubcsat -alg ddfw -cutoff 10000000 -seed $SEED
#./color $GRAPH $SIZE tmp-$$.clq | cadical
