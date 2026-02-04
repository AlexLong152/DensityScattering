#!/bin/bash
make clean
make
# if no arguments then run input.dat
if [ -z "$1" ]; then
  ./onebodyvia1Ndensity input.dat
else
  # Argument provided, use input-<arg>.dat
  ./onebodyvia1Ndensity input-$1.dat
fi
