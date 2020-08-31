#!/bin/bash

source ./setup.sh

BC=`which bc 2> /dev/null`

if [ -z "${BC}" ]; then
  echo "Error: 'bc' not available on your system :("
  exit 1
fi

if [ -z "$1" ]; then
  echo "Usage:"
  echo "  $0 <rawFile_XxYxZ_type.raw>"
  exit 1
fi

IFS='_'
read -rea ADDR <<< "$1"

WORDCOUNT=${#ADDR[@]}

TYPE_POS=`echo "$WORDCOUNT - 1" | $BC`
DIMENSION_POS=`echo "$WORDCOUNT - 2" | ${BC}`
for((I=O; I < $WORDCOUNT - 3; I++)); do
  NAME="${NAME}${ADDR[$I]}_"
done
LAST_WORD_POS=`echo "$WORDCOUNT - 3" | ${BC}`
NAME="${NAME}${ADDR[$LAST_WORD_POS]}"

DIMENSIONS=${ADDR[$DIMENSION_POS]}
TYPE=${ADDR[$TYPE_POS]}

IFS="x"
read -rea ADDR <<< "${DIMENSIONS}"
X=${ADDR[0]}
Y=${ADDR[1]}
Z=${ADDR[2]}

IFS="."
read -rea ADDR <<< "$TYPE"
TYPE=${ADDR[0]}

case $TYPE in 
  "uint8")
    VTKTYPE="unsigned char"
    ;;
  "uint16")
    VTKTYPE="unsigned short"
    ;;
  "int16")
    VTKTYPE="short"
    ;;
  "float32")
    VTKTYPE="float"
    ;;
  "float64")
    VTKTYPE="double"
    ;;
esac 

echo "Converting ${1} into:"
echo "  ${NAME}.vti"
echo "    ${X}x${Y}x${Z}"
echo "    $VTKTYPE"
echo "Processing..."

X=`echo "$X - 1" | ${BC}`
Y=`echo "$Y - 1" | ${BC}`
Z=`echo "$Z - 1" | ${BC}`

$PVPYTHON_DIR/pvpython raw2vtk.py "$1" ${NAME}.vti $X $Y $Z $VTKTYPE
