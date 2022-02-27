GRAPH=$1

./clq.sh $GRAPH > tmp-$$.clq
DSATUR=`./DSATUR/color -filein $GRAPH | tail -n 1 | awk '{print $4}'`

./inner-loop.sh $GRAPH tmp-$$.clq $DSATUR 1
