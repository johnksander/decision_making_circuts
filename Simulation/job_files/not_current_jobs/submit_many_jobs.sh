#!/bin/bash

jobFN=fastD_parsweep.sh

start_idx=1001
total_jobs=50000
jobsz=1000


for idx in $(seq $start_idx $jobsz $total_jobs)
do 	
	sbatch --array=${idx}-$(($idx + $jobsz - 1)) ${jobFN}
	sleep 5
done
