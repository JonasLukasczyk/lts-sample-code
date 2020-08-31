#!/bin/bash

function benchDataSet(){

    THRESHOLD=$1
    LOG_FILE=$2
    OUTPUT=$3
    i=$4

    # Baseline
    echo "    Baseline approach..."
    $PVPYTHON_DIR/pvpython performance_simplification.py $i ${THRESHOLD} 1 0 > ${LOG_FILE}
    DATA_SIZE=`getDataSize ${LOG_FILE}`
    BASELINE_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`
    
    # LTS sequential
    echo "    LTS (1 core)..."
    $PVPYTHON_DIR/pvpython performance_simplification.py $i ${THRESHOLD} 0 1 > ${LOG_FILE}
    LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    # LTS parallel
    echo "    LTS (all cores)..."
    $PVPYTHON_DIR/pvpython performance_simplification.py $i ${THRESHOLD} 1 1 > ${LOG_FILE}
    PARALLEL_LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    echo "$DATA_SIZE $i $BASELINE_AVERAGE_TIME $LTS_AVERAGE_TIME $PARALLEL_LTS_AVERAGE_TIME" >> $OUTPUT

}

function runBench(){

  THRESHOLD=$1
  LOG_FILE="run.log"

  OUTPUT="performance_${THRESHOLD}_simplification.csv"

  echo "Running performance benchmark (persistence threshold: ${THRESHOLD})..."

#  LIST="marmoset_neurons_downsampled.vti random.vti"
  for i in *.vt?; do
#  for i in $LIST; do
    echo "  Processing data set $i..."
    benchDataSet $THRESHOLD $LOG_FILE $OUTPUT $i
  done

  # sort the output file by data size
  cp $OUTPUT output.tmp
  sort -k1 -n output.tmp > $OUTPUT
}

source ./setup.sh

source ./parsing.sh

runBench 1

runBench 99.999
