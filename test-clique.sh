GRAPH=$1
BASE=${GRAPH##*/}
BASE=${BASE%.*}

DIR=tmp
mkdir -p $DIR

echo $GRAPH" "$BASE

NV=`head $GRAPH | awk '/edge/ {print $3}'`
NE=`head $GRAPH | awk '/edge/ {print $4}'`

SAT="SATISFIABLE"
UNS="UNSATISFIABLE"
UNK="UNKNOWN"

CL=`(timeout 1 ./inc_max_clique $GRAPH 1 | tail -n 4 | grep 'Max Clique Size: ' | sed 's/Max Clique Size: //g')`
if [ "$CL" == "" ]; then
  CL=1
fi
echo "c initial clique by cliquer: "$CL

./color $GRAPH $CL | ./cadical/build/cadical | grep "SATIS"
