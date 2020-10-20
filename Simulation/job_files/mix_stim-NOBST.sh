#!/bin/bash
#SBATCH -J mix_stim # A single job name for the array
#SBATCH --time=3-00:00:00 # Running time 
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 750 # Memory request (in Mb)
#SBATCH --account=paul-lab
#SBATCH --partition=paul-compute,neuro-compute,guest-compute
#SBATCH --qos=medium
#SBATCH -o logs/log_parsweep_%A_%a.out
#SBATCH -e logs/log_parsweep_%A_%a.err


cd /work/jksander/rotation/Simulation/
module load share_modules/MATLAB/R2019a

matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "nets_mixstim_NOBSTEST_driver"
