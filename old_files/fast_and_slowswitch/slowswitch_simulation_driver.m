clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'JK';
config_options.sim_name = 'slowswitch_baseline';
options = set_options(config_options);
options.tmax = 2500; %trial simulation time (s) 
stimA = zeros(24,1); %let's do 2 runs per core at 2.5k seconds, try and get 10k switches 
stimB = zeros(size(stimA)); %current modifiers for stimuli (col vectors!)
options.trial_stimuli = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
options.data_recordwindow = [0 999999]; %only record switches in this window (seconds)
options.parforlog = 'on'; 
%---parfor--------------------
num_workers = 12; %parpool workers
addons{1} = fullfile(options.helper_funcdir,'test4switch.m');
addons{2} = fullfile(options.helper_funcdir,'timepoint_stimulus.m');
addons{3} = fullfile(options.helper_funcdir,'next_timepoint.m');
addons{4} = which('poissrnd.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',addons)
%---run-----------------------
modelfile = slowswitch_model(options);
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(gcp('nocreate'))

