clear 
clc
format compact 




%my model 
%---setup---------------------
jID = str2num(getenv('SLURM_ARRAY_TASK_ID'));
t = 5000; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','hpc',...
    'sim_name','parsweep_fastD_baseline','jobID',jID,'tmax',t,...
    'stim_pulse',[t,0],'cut_leave_state',t,'sample_Estay_offset',0);
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)