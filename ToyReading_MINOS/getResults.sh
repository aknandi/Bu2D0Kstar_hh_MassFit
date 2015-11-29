#!/bin/bash

## Loop over individual subjobs
INPUTFILES=()
INPUTDIR="../"
SUBDIR="../"
MIN=1
MAX=10000
GAP=50

for (( i=$MIN; i <= $MAX; i+=$GAP ))
do
  INPUTFILES+=("${INPUTDIR}${SUBDIR}toySeed_$i/toy_$i.txt")
done

## ... Or manually look up results
#INPUTFILES=("../TOYS/test2/toy_31.txt"
#            )

## Read in textfile with floating variables in this toy
readarray VARS < $1
echo ${VARS[@]}

rm results/*.txt

re='^-?[0-9]+([.][0-9]+)?$'

for VAR in ${VARS[*]}
do
  # Only do stuff if not a number (i.e. use the var names only)
  if ! [[ ${VAR} =~ $re ]] ; then
    for FILE in ${INPUTFILES[*]}
    do
      grep ${VAR} ${FILE} >> results/${VAR}.txt
    done
  echo "%s/${VAR} //
  w
  q
  " | ex results/${VAR}.txt
  wc -l results/${VAR}.txt
  fi
done
