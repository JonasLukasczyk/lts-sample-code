#!/bin/bash

# directory where the binary 'pvpython' has been installed
export PVPYTHON_DIR="/usr/local/bin"

# If pvpython fails at detecting TTK, please adjust the following environment
# variable, specifying the directory where the TTK ParaView modules have been 
# installed:
export PV_PLUGIN_PATH="/usr/local/bin/plugins/"


if [ -z "$PVPYTHON_DIR" ]; then
  echo "Error: empty variable PVPYTHON_DIR (installation path for 'pvpython')"
  exit
fi
