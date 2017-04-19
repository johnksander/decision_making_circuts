clear 
clc
format compact 




%paul's model 
%---setup---------------------
config_options.sim_name = 'find_baseline';
JKoptions = set_options(config_options);
JKoptions.current_pulse = 'on';  %'on' | 'off' switch for adding transient pulses (& inital pulse)
JKoptions.tmax = 1000; %trial simulation time (s) 
stimA = [.5:.02:1.5]'; stimB = [.5:.02:1.5]'; %current modifiers for stimuli (col vectors!)
JKoptions.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
JKoptions.parforlog = 'on'; 
%---parfor--------------------
num_workers = 32; %parpool workers
addons = fullfile(JKoptions.helper_funcdir,'test4switch.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{addons})
%---run-----------------------
JK_switching_model(JKoptions)
%---cleanup-------------------
delete(gcp('nocreate'))

