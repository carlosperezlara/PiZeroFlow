#!/bin/bash

#for X in `cat ../runs.dat`
#do
#    rm out/run${X}.root
#    hadd out/run${X}.root out/out_${X}_*root
#    root -b -l -q ComputeTimeConstants.C\($X\);
#done
rm out/all.root
hadd out/all.root out/run*root
