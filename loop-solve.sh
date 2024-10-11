#!/bin/bash
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

#CL=`(timeout 1 ./cliquer/cl $GRAPH -r degree | grep max | tail -n 1 | tr ")" " " | awk '{print $3}')`
CL=`(timeout 1 ./inc_max_clique $GRAPH 1 | tail -n 4 | grep 'Max Clique Size: ' | sed 's/Max Clique Size: //g')`
if [ "$CL" == "" ]; then
  CL=1
fi
echo "c initial clique by cliquer: "$CL

EXT=10
echo -n "c "

for i in $(eval echo "{$CL..300}")
do
  RESULT=`./maxclique $GRAPH $i | ./cadical/build/cadical -c 10000 | grep -e SATIS -e UNKNOWN | awk '{print $2}'`
  if [ "$RESULT" = "$SAT" ]; then
    echo -n "S"
  elif [ "$RESULT" = "$UNS" ]; then
    MAX=$(($i - 1))
    echo "* "$MAX
    ./maxclique $GRAPH $MAX | ./cadical/build/cadical | grep "^v" | ./strip.sh | awk '{if ($1 <= '$NV') print $0}' > $DIR/clique-$$.ord
    cp $DIR/clique-$$.ord cliques/$BASE.clq
    ./optimize-ord $GRAPH $DIR/clique-$$.ord > $DIR/opt-$$.ord
    cp $DIR/opt-$$.ord $DIR/$BASE-$MAX.ord
    rm $DIR/clique-$$.ord
    break
  elif [ "$RESULT" = "$UNK" ]; then
    for j in {1..10}
    do
      echo "" > $DIR/tmp-$$.mod
      echo -n "|"
      UP=$(($i + $j*$EXT))
      CO=$(($j*$EXT*100000))
#      echo $UP" "$CO
#    ./color $GRAPH $UP | yalsat | grep SATIS
      ./color $GRAPH $UP > $DIR/tmp-$$.cnf
#    ~/yalsat-044/palsat tmp-$$.cnf -v $RANDOM
#    ubcsat -i tmp-$$.cnf -alg ddfw -cutoff $CO | grep -v -e "#" -e "=" | awk '{if (NF > 0) print $0}'
      ./ubcsat/ubcsat -i $DIR/tmp-$$.cnf -alg ddfw -cutoff $CO -solve | grep -v -e "#" -e "=" | grep "v " | \
                tr " " "\n" | awk '{if ($1 > 0) print $0}' | grep -v "v" > $DIR/tmp-$UP-$$.mod
      SIZE=`wc $DIR/tmp-$UP-$$.mod | awk '{print $1}'`
#      echo "c SIZE "$SIZE
      if (( "$SIZE" > "1" )); then
        EXT=$(($j*$EXT + 1))
        break
      fi
      rm $DIR/tmp-$$.mod
    done
    for j in $(eval echo "{$EXT..1}")
    do
      SUB=`./maxclique $GRAPH $UP $DIR/tmp-$UP-$$.mod $j | ./cadical/build/cadical | grep SATIS | awk '{print $2}'`
      if [ "$SUB" = "$SAT" ]; then
        echo -n "S"
      elif [ "$SUB" = "$UNS" ]; then
        MAX=$(($UP - $j - 1))
        echo "* "$MAX
        j=$(($j + 1))
        ./maxclique $GRAPH $UP $DIR/tmp-$UP-$$.mod $j | ./cadical/build/cadical | grep "^v" | ./strip.sh | awk '{if ($1 <= '$NV') print $0}' > $DIR/clique-$$.ord
        cp $DIR/clique-$$.ord cliques/$BASE.clq
        ./optimize-ord $GRAPH $DIR/clique-$$.ord > $DIR/opt-$$.ord
        cp $DIR/opt-$$.ord $DIR/$BASE-$MAX.ord
        rm $DIR/clique-$$.ord
        rm $DIR/tmp-$UP-$$.mod
        break
      else
        echo "ERROR"
        break
      fi
    done
    break
  else
    echo "ERROR in clique loop"
    break
    break
  fi
done
#echo $MAX

LEVEL=`./color $GRAPH $MAX $DIR/opt-$$.ord | ./cadical/build/cadical --unsat -c 100000 | grep -e "SATIS" -e "c \- " | tail -n 1 | awk '{print $5}'`
#echo $MAX
#wc opt-$$.ord
if [ "$LEVEL" = "" ]; then
  LEVEL=0
fi
#echo $LEVEL

if (( "$LEVEL" > "15" )); then
  ./color $GRAPH $MAX > $DIR/formula-$$.cnf
  ./ubcsat/ubcsat -i $DIR/formula-$$.cnf -alg ddfw -solve -cutoff 100000000 | grep " 1 1"
  rm $DIR/formula-$$.cnf
  echo "c max clique is equal to chromatic number"
else
  echo -n "c "
  for i in $(eval echo "{$MAX..200}")
  do
#    ./color $GRAPH $i opt-$$.ord | ./cadical/build/cadical --forcephase=1 --phase=1
    RESULT=`./color $GRAPH $i $DIR/opt-$$.ord | ./cadical/build/cadical --unsat --forcephase=1 --phase=0 | grep SATIS | awk '{print $2}'`
    if [ "$RESULT" = "$UNS" ]; then
      echo -n "U"
    elif [ "$RESULT" = "$SAT" ]; then
      echo "* "$i
      if [ "$i" = "$MAX" ]; then
        echo "c max clique is equal to chromatic number"
      fi
      break
   else
      echo "ERROR in bound loop"
      break
    fi
  done
fi

echo -n "c total runtime: "
rm $DIR/opt-$$.ord
