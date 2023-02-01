#!/bin/bash

set -e

g++ -O2 -o A A.cpp
for i in $(seq 0 99)
do
  f=$(printf %04d $i).txt
  echo $f
  ./A < tools/in/$f > out/$f
done
