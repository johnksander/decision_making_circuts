clear 
clc
format compact 




%paul's model 
%---setup---------------------
config_options.sim_name = 'find_baseline';
PMoptions = set_options(config_options);
PMoptions.equal_pools = 'on'; %'on' | 'off'  set switch & stay to have equal properties (b)
PMoptions.current_pulse = 'on';  %'on' | 'off' switch for adding transient pulses (& inital pulse)
PMoptions.tmax = 1000 * 1000; %trial simulation time (ms) 
stimA = [.5:.02:1.5]'; stimB = [.5:.02:1.5]'; %current modifiers for stimuli (col vectors!)
PMoptions.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
PMoptions.parforlog = 'on'; 
%---parfor--------------------
num_workers = 12; %parpool workers
addons = fullfile(PMoptions.helper_funcdir,'test4switch.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{addons})
%---run-----------------------
PM_switching_model(PMoptions)
%---cleanup-------------------
delete(gcp('nocreate'))

