#!/bin/bash

source ./setup.sh

echo "Generating random data set..."

$PVPYTHON_DIR/pvpython random.py
