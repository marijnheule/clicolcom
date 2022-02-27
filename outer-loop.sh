GRAPH=$1

timeout 1 ./inc_max_clique $GRAPH | grep "Clique Output:" | tail -n 1 | sed 's|Clique Output: ||' | tr " " "\n" > tmp-$$.clq
DSATUR=`./DSATUR/color -filein $GRAPH | tail -n 1 | awk '{print $4}'`

./inner-loop.sh $GRAPH tmp-$$.clq $DSATUR 1

rm tmp-$$.clq
