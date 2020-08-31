#!/bin/bash

function benchDataSet(){

    THRESHOLD=$1
    OUTPUT=$2
    i=$3

    LOG_FILE="runRandom_baseline_${THRESHOLD}.log"

    # Baseline
    echo "    Baseline approach..."
    $PVPYTHON_DIR/pvpython performance_randomSimplification.py $i ${THRESHOLD} 1 0 > ${LOG_FILE}
    BASELINE_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`
   
    LOG_FILE="runRandom_lts_${THRESHOLD}.log"
    # LTS sequential
    echo "    LTS (1 core)..."
    $PVPYTHON_DIR/pvpython performance_randomSimplification.py $i ${THRESHOLD} 0 1 > ${LOG_FILE}
    LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    LOG_FILE="runRandom_plts_${THRESHOLD}.log"
    # LTS parallel
    echo "    LTS (all cores)..."
    $PVPYTHON_DIR/pvpython performance_randomSimplification.py $i ${THRESHOLD} 1 1 > ${LOG_FILE}
    PARALLEL_LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    echo "$THRESHOLD $BASELINE_AVERAGE_TIME $LTS_AVERAGE_TIME $PARALLEL_LTS_AVERAGE_TIME" >> $OUTPUT

}

function runBench(){

  THRESHOLD=$1

  OUTPUT="performance_randomSimplification.csv"

  echo "Running performance benchmark (pair threshold: ${THRESHOLD})..."

  INPUT="random.vti"

  benchDataSet $THRESHOLD $OUTPUT $INPUT
}

source ./setup.sh

source ./parsing.sh

for((I=0; I<=100; I+=5)); do
  THRESHOLD=$I
  if [ "$I" == "100" ]; then
    THRESHOLD="99.999"
  fi
  runBench $THRESHOLD
done
