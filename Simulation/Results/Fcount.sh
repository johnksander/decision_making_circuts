#!/usr/bin/env bash

DIR=$1
Ftype=$2 #should be file extension like, mat for *.mat
curr_dir=$(pwd)
Fout=filecount-$DIR
Fout="${Fout//\//-}" #replace directory slashes
(N=$(find $DIR -maxdepth 1 -type f -name '*.'$Ftype | wc -l); cd $curr_dir; echo "number $2 files = $N" >> $Fout) &
