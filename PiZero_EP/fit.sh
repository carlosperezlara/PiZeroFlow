#!/bin/bash

for X in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
    root -b -l -q fit.C\($X,\"$1\"\)
done
