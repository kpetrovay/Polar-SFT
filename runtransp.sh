#!/bin/bash
mkdir res_case1
mkdir res_case2
mkdir res_case3
# 7 x 15 x 11 = 1105 cases/flow type
for flow in 2
do
for u0 in 5 7.5 10 12.5 15 17.5 20
do
for eta in 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750  
do
for tau in 2 3 4 5 6 7 8 9 10 20 100 
do
./transp.py $u0 $eta $tau $flow
done
done
done
done
