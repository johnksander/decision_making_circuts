#!/bin/bash
#SBATCH -J very_slow_sweep # A single job name for the array
#SBATCH --time=3-00:00:00  # Running time 
#SBATCH --mem-per-cpu 750 # Memory request (in Mb)
#SBATCH --partition=ncf
#SBATCH -o /dev/null
#SBATCH -e /dev/null

cd /users/ksander/rotation/Simulation/
module load matlab/R2019a-fasrc01

matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "parsweep_D2t_very_slow_driver"
