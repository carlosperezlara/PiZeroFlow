#!/bin/tcsh

@ N=$1 + 1
set X = `head -$N runs.dat| tail -1`
echo $X

cd BBC_EPC
rm out/run${X}.root
hadd out/run${X}.root out/out_${X}_*root
#root -b -l -q qcent.C\($X\)
#root -b -l -q coef.C\($X\)
#root -b -l -q res.C\($X\)
cd ..
