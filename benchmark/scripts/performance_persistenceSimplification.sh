#!/bin/bash

function benchDataSet(){

    THRESHOLD=$1
    LOG_FILE=$2
    OUTPUT=$3
    i=$4

    # Baseline
    echo "    FTM + Baseline approach..."
    $PVPYTHON_DIR/pvpython performance_persistenceSimplification_ftmBased.py $i ${THRESHOLD}  0 > ${LOG_FILE}
    DATA_SIZE=`getDataSize ${LOG_FILE}`
    FTM_AVERAGE_TIME=`getFTMRunAverageTime ${LOG_FILE}`
    BASELINE_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`
    
    # FTM+LTS
    echo "    FTM + LTS (all core)..."
    $PVPYTHON_DIR/pvpython performance_persistenceSimplification_ftmBased.py $i ${THRESHOLD} 1 > ${LOG_FILE}
    LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    # persistence-LTS
    echo "    persistence-LTS (all cores)..."
    $PVPYTHON_DIR/pvpython performance_persistenceSimplification_ltsBased.py $i ${THRESHOLD} 1 > ${LOG_FILE}
    PERSISTENCE_LTS_AVERAGE_TIME=`getRunAverageTime ${LOG_FILE}`

    echo "$DATA_SIZE $i $FTM_AVERAGE_TIME $BASELINE_AVERAGE_TIME $LTS_AVERAGE_TIME $PERSISTENCE_LTS_AVERAGE_TIME" >> $OUTPUT

}

function runBench(){

  THRESHOLD=$1
  LOG_FILE="runPersistence.log"

  OUTPUT="performance_${THRESHOLD}_persistenceSimplification.csv"

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
