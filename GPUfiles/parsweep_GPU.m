clear 
clc
format compact 


%---setup---------------------
jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 1000; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','hpc',...
    'sim_name','parsweep_D2t_GPU','jobID',jID,'tmax',t,...
    'ratelim_check','on','cut_leave_state',t);
%---run-----------------------
modelfile = spikeout_model_GPU(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)

