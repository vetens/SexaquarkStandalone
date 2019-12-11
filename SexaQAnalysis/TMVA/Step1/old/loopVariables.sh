#!/bin/bash

#dirName = ./LoopCut20Var2016

mkdir ./LoopCut9Var2016
for i in {0..9}
do
   #mkdir ./LoopCut20Var2016/$i
   python BDT_2016_Loop.py $i "LoopCut9Var2016"
done
