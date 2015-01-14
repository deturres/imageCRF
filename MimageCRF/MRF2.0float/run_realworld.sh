#!/bin/bash

# script to run multiple images and multiple values 
# of a and b

echo "Hello World"

neg1=`expr -1`

for a in 1 2 5 10
do
  for b in -10 -5 -2 -1 1 2 5 10
  do
    echo "a = $a, b = $b"
    for file in /localdisk/videosegmentation/NIPSexp/noise/*
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
      fileout="/localdisk/videosegmentation/NIPSexp/noiseout/$bname2$astr$bstr"
      echo $fileout
      ./example_dd_realworld $file $fileout $a $b
    done   
  done
done