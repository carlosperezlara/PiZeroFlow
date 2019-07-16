#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segmentsERT.dat| tail -1`
echo $A
@ B=-1

./Run_PiZero_EP $A $B ERT
./Run_PiZero_EP $A $B ERTD0
./Run_PiZero_EP $A $B ERTD1
./Run_PiZero_EP $A $B ERTA0
./Run_PiZero_EP $A $B ERTA1
./Run_PiZero_EP $A $B ERTT0
./Run_PiZero_EP $A $B ERTT1
./Run_PiZero_EP $A $B ERTFD0
./Run_PiZero_EP $A $B ERTFD1
./Run_PiZero_EP $A $B ERTFA0
./Run_PiZero_EP $A $B ERTFA1
./Run_PiZero_EP $A $B ERTFT0
./Run_PiZero_EP $A $B ERTFT1


exit
