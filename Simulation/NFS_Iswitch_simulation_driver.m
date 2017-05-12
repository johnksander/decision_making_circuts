clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'JK';
config_options.sim_name = 'NFSv2_Iswitch';
options = set_options(config_options);
options.current_pulse = 'off';  %'on' | 'off' set off for low memory implementation
options.tmax = 5000; %trial simulation time (s) 
stimA = ones(24,1); %let's do 1 run per core at 10k seconds, try and get 10k switches 
stimA = stimA + .5; %set stim A bias to 1.5
stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
%stimA = [1:.02:1.5]'; stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
options.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
options.NFS = 'Iswitch'; %cells driving forced switch:  Estay | Eswitch | Istay | Iswitch
options.NFS_onset_min = 1; %minimum state duration for forced switch (seconds)
options.NFS_stoppush = 200e-3; %noise push duration (ms) 
options.NFS_noisepush = 100e-12; %push noise by 100 picoamps 
options.parforlog = 'on'; 
%---parfor--------------------
num_workers = 12; %parpool workers
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
