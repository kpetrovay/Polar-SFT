#!/bin/sh
mkdir maps_case1
mkdir maps_case2
mkdir maps_case3
for flow in 2
do
for tau in 2 3 4 5 6 7 8 9 10 20 100 
do
#echo "Now running check.py with args" `expr $u0` `expr $eta` `expr $tau` `expr $flow`
./parammap.py $flow $tau
done
done
