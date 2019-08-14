#!/bin/tcsh

@ N=$1 + 1
set X = `head -$N runs.dat| tail -1`
echo $X

#cd BBC_EPC
#rm out/run${X}.root
#hadd out/run${X}.root out/out_${X}_*root
#root -b -l -q qcent.C\($X\)
#root -b -l -q coef.C\($X\)
#root -b -l -q res.C\($X\)
#cd ..

#cd Basic
#rm out/run${X}.root
#hadd out/run${X}.root out/out_${X}_*root
#cd ..

#cd BBC_EPC
#rm out/run${X}.root
#hadd out/run${X}.root out/out_${X}_*root
#cd ..

cd PiZero_EP
foreach Y ( "A0" "A1" "D0" "D1" "FA0" "FA1" "FD0" "FD1" "FT0" "FT1" "T0" "T1" )
#foreach Y ( "NOM" )
    rm out${Y}/run${X}.root
    hadd out${Y}/run${X}.root out${Y}/out_${X}*.root
end
cd ..

#cd EventChecker
#rm out/run${X}.root
#hadd out/run${X}.root out/out_${X}_*root
#root -b -l -q ComputeTimeConstants.C\($X\);
#cd ..
