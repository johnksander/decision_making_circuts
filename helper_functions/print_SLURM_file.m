function print_SLURM_file(options,job_runtime,partition)
%creates a SLURM batch file 
%job_runtime is in hours 

%partition = 'paul'; %guest or paul

get_err = 'no'; %'yes' | 'no' ,for printing error output 

M = mod(job_runtime,1);
M = round(M * 60);
H = floor(job_runtime);

walltime = sprintf('%02.0f:%02.0f:00',H,M);

outFN = fullfile(options.batchdir,sprintf('batch_job_%s.sh',partition));

if exist(outFN,'file') > 0,delete(outFN);end
f = fopen(outFN,'w');
fprintf(f,'#!/bin/bash\n'); 
fprintf(f,'#SBATCH -J EQstim_batch\n'); 
fprintf(f,'#SBATCH --time=%s\n',walltime); 
fprintf(f,'#SBATCH --cpus-per-task 1\n');
fprintf(f,'#SBATCH --mem-per-cpu 750\n');
fprintf(f,'#SBATCH --account=paul-lab\n');
switch partition
    case 'paul'
        fprintf(f,'#SBATCH --partition=paul-compute,neuro-compute,guest-compute\n');
        fprintf(f,'#SBATCH --qos=medium\n');
    case 'guest'
        fprintf(f,'#SBATCH --partition=guest-compute\n');
        fprintf(f,'#SBATCH --qos=low\n');
end
fprintf(f,'#SBATCH -o /dev/null\n');
switch get_err
    case 'yes'
        fprintf(f,'#SBATCH -e %s\n',fullfile(options.batchdir,'errout_%%A_%%a.err')); %for error output
    case 'no'
        fprintf(f,'#SBATCH -e /dev/null\n');
end
%okay, now for the bash stuff
fprintf(f,'\n\n\n');


fprintf(f,'cd %s\n',options.batchdir);
fprintf(f,'module load share_modules/MATLAB/R2019a\n\n');

fprintf(f,'matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "driver"\n');
%fprintf(f,'srun matlab -singleCompThread -nodisplay -nodesktop -nosplash -r "driver"\n');


fclose(f);
