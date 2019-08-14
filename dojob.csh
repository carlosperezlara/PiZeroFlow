#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
#set A = `head -$N sss.dat| tail -1`
echo $A
@ B=-1

#./Run_ReadTree $A $B

#./Run_MX_EPC $A $B
#./Run_PiZero $A $B

#./Run_EventChecker $A $B
#./Run_BBC_EPC $A $B

#./Run_PiZero_EP $A $B 
./Run_PiZero_EP $A $B D0
./Run_PiZero_EP $A $B D1
./Run_PiZero_EP $A $B A0
./Run_PiZero_EP $A $B A1
./Run_PiZero_EP $A $B T0
./Run_PiZero_EP $A $B T1
./Run_PiZero_EP $A $B FD0
./Run_PiZero_EP $A $B FD1
./Run_PiZero_EP $A $B FA0
./Run_PiZero_EP $A $B FA1
./Run_PiZero_EP $A $B FT0
./Run_PiZero_EP $A $B FT1

#./Run_Charged_QC $A $B

exit
