clear 
clc
format compact 




%my model 
%---setup---------------------
%jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
jID = str2num(getenv('SLURM_ARRAY_TASK_ID'));
t = 5000; %trial simulation time (s) 
options = set_options('modeltype','NETS','comp_location','hpc',...
    'timestep',.02e-3,...
    'sim_name','dttest_0-02_fastD','jobID',jID,'tmax',t,...
    'stim_pulse',[t,0],'cut_leave_state',t,'sample_Estay_offset',0);
%---run-----------------------
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