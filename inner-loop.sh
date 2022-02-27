#set -x
GRAPH=$1  # DIMACS file
LOWER=$2  # clique file
UPPER=$3  # number
DEPTH=$4  # number

SAT="SATISFIABLE"
UNS="UNSATISFIABLE"
UNK="UNKNOWN"

CUTOFF=$(($DEPTH*1000000))

cp $LOWER tmp-$$.clq

CLIQUE=`wc $LOWER | awk '{print $1}'`
echo "c initial clique has size "$CLIQUE
MAX=$UPPER
for BOUND in $(eval echo "{$MAX..0}")
do
  echo $BOUND
  ./color $GRAPH $BOUND $LOWER | ./ubcsat/ubcsat -alg walksat -solve -cutoff $CUTOFF > walk-$$.txt
  if grep  " 1 1 " walk-$$.txt
  then
    MAX=$BOUND
    cat walk-$$.txt | ./strip.sh > tmp-$$.col
  else
    break
  fi
done
echo $MAX

echo -n "c "
OFFSET=$(($MAX-$CLIQUE-1))
for i in $(eval echo "{$OFFSET..0}")
do
  ./maxclique $GRAPH $MAX tmp-$$.col $i | ./cadical/build/cadical --unsat --forcephase=1 --phase=0 -c $CUTOFF > cdcl-$$.tmp
  RESULT=`cat cdcl-$$.tmp | grep SATIS | awk '{print $2}'`
  if [ "$RESULT" = "$SAT" ]; then
    NEWBOUND=$(($MAX - $i))
    cat cdcl-$$.tmp | ./strip.sh | head -n $NEWBOUND > tmp-$$.clq
    echo -n "S"
  elif [ "$RESULT" = "$UNS" ]; then
    i=$(($i+1))
    SIZE=$(($MAX-$i))
    echo " found clique number: "$SIZE

    ./optimize-ord $GRAPH tmp-$$.clq > tmp-$$.ord
    echo -n "c "
    for j in $(eval echo "{$SIZE..200}")
    do
      RESULT=`./color $GRAPH $j tmp-$$.ord | ./cadical/build/cadical --unsat --forcephase=1 --phase=0 | grep SATIS | awk '{print $2}'`
      if [ "$RESULT" = "$UNS" ]; then
        echo -n "U"
      elif [ "$RESULT" = "$SAT" ]; then
      rm tmp-$$.ord
      echo " found chromatic number: "$j
      break
    fi
    done

    break
  else
    echo
    DEPTH=$(($DEPTH*2))
    ./inner-loop.sh $GRAPH tmp-$$.clq $MAX $DEPTH
    break
  fi
done
