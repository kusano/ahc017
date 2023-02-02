#!/bin/bash

set -e

g++ -O2 -o A -DTOPCODER_LOCAL=1 A.cpp
for i in $(seq 0 99)
do
  f=$(printf %04d $i).txt
  echo -n $f
  ./A < tools/in/$f > out/$f
done
