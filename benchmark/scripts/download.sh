#!/bin/bash

WGET=`which wget 2> /dev/null`

if [ -z "$WGET" ]; then
  echo "Error: 'wget' not available on your system :("
  exit
fi

for i in `cat url.txt`; do 
  $WGET $i
done
