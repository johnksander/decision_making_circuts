clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'PS_stim';
config_options.sim_name = 'parsweep_stims_run2';
config_options.jobID = str2num(getenv('SGE_TASK_ID'));
config_options.tmax = 5000; %trial simulation time (s) 
config_options.force_back2stay = true;
options = set_options(config_options);
%---stimulus & network---------

%this is all configured in set_options

%options.stim_targs = 'baseline'; % Eswitch | Estay | baseline
%Rstim = 0; %stimulus (in Hz) to target cells 
%stimA = Rstim; %one trial per job 
%stimB = Rstim; %one trial per job 
%give the same to stim B, just so it doesn't get stuck forever 
%options.trial_stimuli = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
%---run-----------------------
modelfile = parsweep_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
