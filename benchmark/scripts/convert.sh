#!/bin/bash

source ./setup.sh

for i in *raw; do 
  ./raw2vtk.sh $i
done
