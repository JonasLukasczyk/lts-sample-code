#!/bin/bash

source ./setup.sh

if [ ! -e oceanVortices.vtu ]; then
  cp oceanVortices.vtu.old oceanVortices.vtu
fi

echo "Pre-processing the 'oceanVortices.vtu' data set..."

$PVPYTHON_DIR/pvpython preProcessOceanVortices.py

mv oceanVortices.vtu oceanVortices.vtu.old
