#!/bin/bash
#SBATCH -J optimStim # A single job name for the array
#SBATCH --time=3-00:00:00  # Running time 
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 1000 # Memory request (in Mb)
#SBATCH --account=paul-lab
#SBATCH --partition=paul-compute,neuro-compute,guest-compute
#SBATCH --qos=medium
#SBATCH -o log_optimStim_%A_%a.out  #or log_parsweep_%A_%a.out
#SBATCH -e log_optimStim_%A_%a.err
#SBATCH --array=1-1 #without an array ID the random number seed gen fails

cd /work/jksander/rotation/Simulation/
module load share_modules/MATLAB/R2019a

matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "equate_D2t_driver"
