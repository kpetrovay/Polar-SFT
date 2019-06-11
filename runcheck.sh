#!/bin/sh
mkdir plots_case1
mkdir plots_case2
mkdir plots_case3
for flow in 2 
do
rm -f params`expr $flow`.dat
for u0 in 5 7.5 10 12.5 15 17.5 20
do
for eta in 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750  
do
for tau in 2 3 4 5 6 7 8 9 10 20 100 
do
echo "Now running check.py with args" `expr $u0` `expr $eta` `expr $tau` `expr $flow`
./check.py $u0 $eta $tau $flow
done
done
done
done
