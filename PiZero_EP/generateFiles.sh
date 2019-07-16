#!/bin/bash

#for Y in NOM A0  A1  D0  D1  FA0 FA1 FD0 FD1 FT0 FT1 T0  T1
for Y in NOM
do
    rm allfiles/all${Y}.root
    hadd allfiles/all${Y}.root  out${Y}/run*root
done

