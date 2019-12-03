#!/bin/bash
#SBATCH -J D2t # A single job name for the array
#SBATCH --time=24:00:00	 # Running time 
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 1750 # Memory request (in Mb)
#SBATCH --account=paul-lab
#SBATCH --partition=guest-compute
#SBATCH --qos=low
#SBATCH -o /dev/null  #or log_parsweep_%A_%a.out
#SBATCH -e /dev/null


cd /work/jksander/rotation/Simulation/
module load share_modules/MATLAB/R2019a

matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "parsweep_D2t_driver"
