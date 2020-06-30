#!/usr/bin/env bash

JFILE=$1
#JFILE=driver_ItoEfixed_2
#JFILE=driver_EtoIfixed_3
LOG=job_out_$JFILE.txt

module load MATLAB/R2019a
nohup matlab -nodisplay -nodesktop -nosplash -r "$JFILE;exit" >> $LOG 2>&1 &
