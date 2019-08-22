#!/bin/bash

for Y in NOM A0 A1 D0 D1 T0 T1 FA0 FA1 FD0 FD1 FT0 FT1
do
    for X in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
    #for X in 10 11 12 13 14 15 16 17 18
    #for X in 15 16 17 18
    #for X in 18
    do
	root -b -l -q fit.C\($X,\"EP\",\"$Y\"\)
    done
done
