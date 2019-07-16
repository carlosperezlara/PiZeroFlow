#!bash

rm runs.info.dat
for X in `cat runs.dat`
do
    ./getDBinfo $X >> runs.info.dat
done
