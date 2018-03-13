clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'PS';
config_options.sim_name = 'parsweep_baseline';
config_options.jobID = str2num(getenv('SGE_TASK_ID'));
%config_options.jobID = 5;
config_options.tmax = 5000; %trial simulation time (s) 
options = set_options(config_options);
%---stimulus------------------
options.stim_targs = 'baseline'; % Eswitch | Estay | baseline
Rstim = 0; %stimulus (in Hz) to target cells 
stimA = Rstim; %one trial per job 
stimB = Rstim; %one trial per job 
%give the same to stim B, just so it doesn't get stuck forever 
options.trial_stimuli = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
%---run-----------------------
modelfile = parsweep_model(options);
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)

