#!/bin/bash

for X in `cat ../runs.dat`
do
    rm out/run${X}.root
    hadd out/run${X}.root out/out_${X}_*root
    #root -b -l -q qcent.C\($X\)
    #root -b -l -q coef.C\($X\)
    #root -b -l -q res.C\($X\)
done
