#!/bin/bash
#SBATCH -J job_name # A single job name for the array
#SBATCH --time=24:00:00	 # Running time 
#SBATCH --account=paul-lab
#SBATCH --partition=paul-gpu
#SBATCH --qos=medium
#SBATCH --gres=gpu:GTX:1   #do I have to specify the kind of GPU, or can I just take any?
#SBATCH --nodes=1 #one node per task
#SBATCH --array=1-100 #many jobs 

cd /work/jksander/rotation/Simulation/
module load share_modules/MATLAB/R2019a

matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "parsweep_GPU"


