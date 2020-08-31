#!/bin/bash

# config
# ==============================================================================

INPUT1="performance_1_simplification.csv"
INPUT2="performance_99.999_simplification.csv"

TIME_DIGITS="3"
SPEEDUP_DIGITS="1"
EFFICIENCY_DIGITS="0"
ALT_COLOR="e1e1e1"

OUTPUT="performance_simplification.tex"

# deps
# ==============================================================================
BC=`which bc 2> /dev/null`

if [ -z "${BC}" ]; then
  echo "Error: 'bc' not available on your system :("
  exit 1
fi

CAT=`which cat 2> /dev/null`

if [ -z "${CAT}" ]; then
  echo "Error: 'cat' not available on your system :("
  exit 1
fi


HEAD=`which head 2> /dev/null`

if [ -z "${HEAD}" ]; then
  echo "Error: 'head' not available on your system :("
  exit 1
fi

NPROC=`which nproc 2> /dev/null`

if [ -z "${NPROC}" ]; then
  echo "Error: 'nproc' not available on your system :("
  exit
fi

PRINTF=`which printf 2> /dev/null`

if [ -z "${PRINTF}" ]; then
  echo "Error: 'printf' not available on your system :("
  exit 1
fi

SED=`which sed 2> /dev/null`

if [ -z "${SED}" ]; then
  echo "Error: 'sed' not available on your system :("
  exit 1
fi

TAIL=`which tail 2> /dev/null`

if [ -z "${TAIL}" ]; then
  echo "Error: 'tail' not available on your system :("
  exit 1
fi


WC=`which wc 2> /dev/null`

if [ -z "${WC}" ]; then
  echo "Error: 'wc' not available on your system :("
  exit 1
fi

# script
# ==============================================================================

# assuming hyper-threading
CORE_NUMBER=`echo $(${NPROC}) / 2.0 | ${BC}`

$SED "s/xx/${CORE_NUMBER}/g" performance_simplification_latexHeader.tex > $OUTPUT

DATASET_NUMBER=`${WC} -l < $INPUT1`

for((I=1; I <= ${DATASET_NUMBER}; I++)); do
  
  # parse line from file 1
  LINE1=`$HEAD -n$I $INPUT1 | $TAIL -n1`
  IFS=" "
  read -rea ADDR <<< ${LINE1}
  DATA_SIZE1=`${PRINTF} "%'.f\n" ${ADDR[0]}`
  DATA_NAME1=${ADDR[1]/.vt?/}
  DATA_NAME1=`tr '[:lower:]' '[:upper:]'  <<< "${DATA_NAME1:0:1}"`${DATA_NAME1:1}
  BASELINE_TIME1=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[2]}`
  LTS_TIME1=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[3]}`
  LTS_SPEEDUP1=`echo "$BASELINE_TIME1/$LTS_TIME1" | $BC -l`
  LTS_SPEEDUP1=`${PRINTF} "%.${SPEEDUP_DIGITS}f" ${LTS_SPEEDUP1}`
  PLTS_TIME1=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[4]}`
  PLTS_SPEEDUP1=`echo "$BASELINE_TIME1/$PLTS_TIME1" | $BC -l`
  PLTS_SPEEDUP1=`${PRINTF} "%.${SPEEDUP_DIGITS}f" ${PLTS_SPEEDUP1}`
  PLTS_EFFICIENCY1=`echo "100*(($LTS_TIME1/$PLTS_TIME1)/$CORE_NUMBER)" | $BC -l`
  PLTS_EFFICIENCY1=`${PRINTF} "%.${EFFICIENCY_DIGITS}f" ${PLTS_EFFICIENCY1}`

  echo "      ${DATA_NAME1/_/\\_} &" >> $OUTPUT
  echo "      $DATA_SIZE1 &" >> $OUTPUT
  echo "      $BASELINE_TIME1 &" >> $OUTPUT
  echo "      $LTS_TIME1 &" >> $OUTPUT
  echo "      \\textbf{$LTS_SPEEDUP1} &" >> $OUTPUT
  echo "      $PLTS_TIME1 &" >> $OUTPUT
  echo "      \\textbf{$PLTS_SPEEDUP1} &" >> $OUTPUT
  echo "      \\emph{$PLTS_EFFICIENCY1 \\%}\\\\" >> $OUTPUT
  if [ "$I" != "${DATASET_NUMBER}" ]; then
    echo "      \\arrayrulecolor{lightgray}\\hline"  >> $OUTPUT
  fi
done

  echo "      \\arrayrulecolor{black}\\hline"  >> $OUTPUT

for((I=1; I <= ${DATASET_NUMBER}; I++)); do
  
  # parse line from file 2
  LINE2=`$HEAD -n$I $INPUT2 | $TAIL -n1`
  IFS=" "
  read -rea ADDR <<< ${LINE2}
  DATA_SIZE2=`${PRINTF} "%'.f\n" ${ADDR[0]}`
  DATA_NAME2=${ADDR[1]/.vt?/}
  DATA_NAME2=`tr '[:lower:]' '[:upper:]'  <<< "${DATA_NAME2:0:1}"`${DATA_NAME2:1}
  BASELINE_TIME2=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[2]}`
  LTS_TIME2=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[3]}`
  LTS_SPEEDUP2=`echo "$BASELINE_TIME2/$LTS_TIME2" | $BC -l`
  LTS_SPEEDUP2=`${PRINTF} "%.${SPEEDUP_DIGITS}f" ${LTS_SPEEDUP2}`
  PLTS_TIME2=`${PRINTF} "%.${TIME_DIGITS}f" ${ADDR[4]}`
  PLTS_SPEEDUP2=`echo "$BASELINE_TIME2/$PLTS_TIME2" | $BC -l`
  PLTS_SPEEDUP2=`${PRINTF} "%.${SPEEDUP_DIGITS}f" ${PLTS_SPEEDUP2}`
  PLTS_EFFICIENCY2=`echo "100*(($LTS_TIME2/$PLTS_TIME2)/$CORE_NUMBER)" | $BC -l`
  PLTS_EFFICIENCY2=`${PRINTF} "%.${EFFICIENCY_DIGITS}f" ${PLTS_EFFICIENCY2}`

  echo "      \\cellcolor[HTML]{${ALT_COLOR}}${DATA_NAME2/_/\\_} &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}}$DATA_SIZE2 &" >> $OUTPUT

  echo "      \\cellcolor[HTML]{${ALT_COLOR}} $BASELINE_TIME2 &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}} $LTS_TIME2 &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}} \\textbf{$LTS_SPEEDUP2} &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}} $PLTS_TIME2 &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}} \\textbf{$PLTS_SPEEDUP2} &" >> $OUTPUT
  echo "      \\cellcolor[HTML]{${ALT_COLOR}} \\emph{$PLTS_EFFICIENCY2 \\%}\\\\" >> $OUTPUT
  echo "      \\arrayrulecolor{lightgray}\\hline"  >> $OUTPUT

done

$CAT performance_simplification_latexFooter.tex >> $OUTPUT
