clear 
clc
format compact 


setenv('SLURM_JOBID','1010')
setenv('SLURM_ARRAY_TASK_ID','999')
%my model 
%---setup---------------------
jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 600; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','woodstock',...
    'sim_name','test','jobID',jID,'tmax',t,...
    'ratelim_check','on','cut_leave_state',t);
%---run-----------------------
 options.ItoE = 12.4541; options.EtoI = 0.2478;
modelfile = spikeout_model(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)
