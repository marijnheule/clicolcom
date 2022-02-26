GRAPH=$1

CDCLLIM=10000000
DDFWLIM=10000000

SEED=123456

SAT="SATISFIABLE"
UNS="UNSATISFIABLE"
UNK="UNKNOWN"

echo $GRAPH
timeout 1 ./inc_max_clique $GRAPH > tmp-clq-$$.txt
cat tmp-clq-$$.txt | grep "Clique Output:" | tail -n 1 | sed 's|Clique Output: ||' | tr " " "\n" > tmp-$$.clq
SIZE=`wc tmp-$$.clq | awk '{print $1}'`
if grep -q terminated tmp-clq-$$.txt
then
  rm tmp-clq-$$.txt
  echo "c found max clique of size "$SIZE
  CUT=$(($CDCLLIM/$SIZE))
  ./color $GRAPH $SIZE tmp-$$.clq | ./cadical/build/cadical --forcephase=1 --phase=0 -c $CUT | grep SATIS > tmp-cdcl-$$.txt
  if grep -q "s SATISFIABLE" tmp-cdcl-$$.txt
  then
    rm tmp-cdcl-$$.txt
    echo "c found chromatic number: "$SIZE
  elif grep -q "s UNSATISFIABLE" tmp-cdcl-$$.txt
  then
    rm tmp-cdcl-$$.txt
    echo "c chromatic number larger than clique number"
    ./optimize-ord $GRAPH tmp-$$.clq > tmp-$$.ord
    SIZE=$(($SIZE+1))
    echo -n "c "
    for i in $(eval echo "{$SIZE..200}")
    do
      RESULT=`./color $GRAPH $i tmp-$$.ord | ./cadical/build/cadical --unsat --forcephase=1 --phase=0 | grep SATIS | awk '{print $2}'`
      if [ "$RESULT" = "$UNS" ]; then
        echo -n "U"
      elif [ "$RESULT" = "$SAT" ]; then
        rm tmp-$$.ord
        echo " found clique number: "$i
        break
      fi
    done
  else
    rm tmp-cdcl-$$.txt
    echo "c timeout, switching to local search"
    ./color $GRAPH $SIZE tmp-$$.clq | ./ubcsat/ubcsat -alg ddfw -cutoff $DDFWLIM -seed $SEED | grep -v "#" | grep " 1 " > tmp-ddfw-$$.txt
    cat tmp-ddfw-$$.txt
    if grep -q " 1 1 " tmp-ddfw-$$.txt
    then
      echo " c found clique number: "$SIZE
    fi
    rm tmp-ddfw-$$.txt
  fi
else
  rm tmp-clq-$$.txt
  echo "c terminated with clique of size "$SIZE
fi

rm tmp-$$.clq
