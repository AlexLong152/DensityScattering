#!/bin/bash
rm *.o *.mod
gfortran -c parseFile.f
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to compile parseFile.f"
  exit 1
fi

gfortran -c PionScat.f
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to compile PionScat.f"
  exit 1
fi

gfortran -c CalcCrossSec.f
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to compile main.f"
  exit 1
fi

gfortran *.o -o PionScat
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to link executable"
  exit 1
fi

./PionScat
