#!/usr/bin/env bash

DIR=$1
#DIR=parsweep_fastD_baseline

nohup zip -r $DIR.zip $DIR >> out4zip_job.txt 2>&1 &
