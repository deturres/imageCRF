#!/bin/bash

# script to run multiple images and multiple values 
# of a and b

echo "Hello World"

neg1=`expr -1`

c=0.5
d=-1

for a in 1
do
  for b in -15 # 0 # 0.5 0.2 -0.2 -0.5 #1 -1 -1.1 -2 -3 -5 -7 -10
  do
    echo "a = $a, b = $b"
    for file in /localdisk/videosegmentation/NIPSexp/FlipperGraphs/ims5/*
    do
      bname=`basename $file`
      bname2=`echo $bname | sed 's/\(.*\)\..*/\1/'`
      # echo $bname2
      astr=a;
      bstr=b;
      if [ $a -lt 0 ]
      then      
        na=`expr $neg1 \* $a`;
        astr="n$na"
      else
        astr="p$a"
      fi
      if [ $b -lt 0 ]
      then      
        nb=`expr $neg1 \* $b`;
        bstr="n$nb"
      else
        bstr="p$b"  
      fi
      fileout="/localdisk/videosegmentation/NIPSexp/FlipperGraphs/ims5out/$bname2$astr$bstr"
      fileouttxt="/localdisk/videosegmentation/NIPSexp/FlipperGraphs/ims5out/$bname2$astr$bstr.txt"
      echo $fileout
      ./example_dd_realworld2 $file $fileout $a $b $c $d > $fileouttxt
    done   
  done
done