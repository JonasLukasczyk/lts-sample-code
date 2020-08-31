#!/bin/bash

BC=`which bc 2> /dev/null`

if [ -z "${BC}" ]; then
  echo "Error: 'bc' not available on your system :("
  exit 1
fi


function getDataSize(){

  WORD_COUNT="0"
  for WORD in `cat $1 | grep "Size of input data"`; do
    WORD_COUNT=`echo "$WORD_COUNT + 1" | ${BC}`
    #echo "Word #${WORD_COUNT}: ${WORD}"  

    # going for the word of interest
    if [ "$WORD_COUNT" == "5" ]; then
      SIZE=${WORD}
    fi
  done

  echo $SIZE
}

function computeAverageRunTime(){

  TIME_ARRAY=$1

  # sort the initial array
  IFS=$'\n' SORTED_TIME_ARRAY=($(sort -n <<<"${TIME_ARRAY[*]}"))

  TIME_SUM="0"
  for((TIME=1; TIME<${#SORTED_TIME_ARRAY[@]} - 1;TIME++)); do
    TIME_SUM=`echo "${TIME_SUM} + ${SORTED_TIME_ARRAY[${TIME}]}" | ${BC} `
    #echo "TIME: ${TIME}: ${SORTED_TIME_ARRAY[${TIME}]} SUM: ${TIME_SUM}"
  done

  AVERAGE_TIME=`echo "${TIME_SUM}/10.0" | ${BC} -l`

  echo "$AVERAGE_TIME"

}

function getLTSRunAverageTime(){

  WORD_COUNT="0"
  TIME_ARRAY=""
  for WORD in `cat $1 | grep "Complete"`; do
    WORD_COUNT=`echo "$WORD_COUNT + 1" | ${BC}`
    
    # going for the line of interest
    if [ $WORD_COUNT == "4" ]; then
      WORD_COUNT="0"
      IFS="|"
      read -rea ADDR <<< ${WORD}
      TIME=${ADDR[0]}
      TIME=${TIME/[/}
      TIME=${TIME/s/}
      TIME_ARRAY=( "${TIME_ARRAY[@]}" "${TIME}" )
    fi
  done

  computeAverageRunTime ${TIME_ARRAY[@]}
}

function getBaseLineRunAverageTime(){
  
  WORD_COUNT="0"
  for WORD in `cat $1 | grep "simplified"`; do
    WORD_COUNT=`echo "$WORD_COUNT + 1" | ${BC}`
    #echo "Word #${WORD_COUNT}: ${WORD}"  

    # going for the word of interest
    if [ "$WORD_COUNT" == "6" ]; then
      TIME=${WORD}
      TIME_ARRAY=( "${TIME_ARRAY[@]}" "${TIME}" )
    fi

    if [ "$WORD_COUNT" == "11" ]; then
      WORD_COUNT="0"
    fi 
  done

  computeAverageRunTime ${TIME_ARRAY[@]}
}

function getFTMRunAverageTime(){
  
  WORD_COUNT="0"
  for WORD in `cat $1 | grep "Total  in"`; do
    WORD_COUNT=`echo "$WORD_COUNT + 1" | ${BC}`
    #echo "Word #${WORD_COUNT}: ${WORD}"  

    # going for the word of interest
    if [ "$WORD_COUNT" == "4" ]; then
      TIME=${WORD}
      TIME_ARRAY=( "${TIME_ARRAY[@]}" "${TIME}" )
      WORD_COUNT="0"
    fi
  done

  computeAverageRunTime ${TIME_ARRAY[@]}
}

function getRunAverageTime(){
  if [ -z "`cat $1 | grep "Complete"`" ]; then
    getBaseLineRunAverageTime $1
  else 
    getLTSRunAverageTime $1
  fi
}
