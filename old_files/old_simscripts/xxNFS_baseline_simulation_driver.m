clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'JK';
config_options.sim_name = 'NFSv5_baseline';
options = set_options(config_options);
options.current_pulse = 'off';  %'on' | 'off' set off for low memory implementation
options.tmax = 5000; %trial simulation time (s) 
stimA = ones(32,1); %let's do 1 run per core at 10k seconds, try and get 10k switches 
stimA = stimA + .1; %set stim A bias to 1.1
stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
%stimA = [1:.02:1.5]'; stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
options.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
options.NFS = 'Eswitch'; %cells driving forced switch:  Estay | Eswitch | Istay | Iswitch
options.NFS_onset_min = .75; %minimum state duration for forced switch (seconds)
options.NFS_recordwindow = [.75 99999]; %only record switches in this window (seconds)
options.NFS_stoppush = max(options.NFS_recordwindow)-options.NFS_onset_min; %noise push duration (s)
options.NFS_noisepush = 0; %dummy noise push for baseline 
options.parforlog = 'on'; 
%---parfor--------------------
num_workers = 16; %parpool workers
addons{1} = fullfile(options.helper_funcdir,'test4switch.m');
addons{2} = fullfile(options.helper_funcdir,'timepoint_current.m');
addons{3} = fullfile(options.helper_funcdir,'next_timepoint.m');
addons{4} = fullfile(options.helper_funcdir,'NFSnoise_adjustment.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',addons)
%---run-----------------------
modelfile = noisycurrent_model(options);
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(gcp('nocreate'))
