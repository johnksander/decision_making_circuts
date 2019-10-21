#!/usr/bin/env bash

export SIM_NAME=$1
export JID=$2

JFILE=inspect
LOG=job_out_$JFILE.txt

module load MATLAB/R2019a
nohup matlab -nodisplay -nodesktop -nosplash -r "$JFILE;exit" >> $LOG 2>&1 &
