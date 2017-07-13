clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'JK';
config_options.sim_name = 'fastswitch_stimulus';
options = set_options(config_options);
options.tmax = 2500; %trial simulation time (s) 
Rstim = 207.5; %207.5 hz stimulus to Estay cells 
stimA = zeros(32,1) + Rstim; %let's do 2 runs per core at 2.5k seconds, try and get 10k switches 
stimB = zeros(size(stimA)) + Rstim; %these must be col vectors!
%give the same to stim B, just so we get 2x as much data 
options.trial_stimuli = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
options.data_recordwindow = [0 999999]; %only record switches in this window (seconds)
options.parforlog = 'on'; 
%---parfor--------------------
num_workers = 16; %parpool workers
addons{1} = fullfile(options.helper_funcdir,'test4switch.m');
addons{2} = fullfile(options.helper_funcdir,'timepoint_stimulus.m');
addons{3} = fullfile(options.helper_funcdir,'next_timepoint.m');
addons{4} = which('poissrnd.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',addons)
%---run-----------------------
modelfile = fastswitch_model(options);
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(gcp('nocreate'))

