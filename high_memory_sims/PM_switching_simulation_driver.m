clear 
clc
format compact 




%Paul's model 
%---setup---------------------
config_options.modeltype = 'PM';
config_options.sim_name = 'switching_bias2';
PMoptions = set_options(config_options);
PMoptions.equal_pools = 'off'; %'on' | 'off' set switch & stay to have equal properties (b)
PMoptions.current_pulse = 'off';  %'on' | 'off' set off for low memory implementation
PMoptions.tmax = 1000 * 1000; %trial simulation time (ms) 
stimA = [1:.02:2]'; stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
PMoptions.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
PMoptions.parforlog = 'on'; 
%---parfor--------------------
num_workers = 12; %parpool workers
addons{1} = fullfile(PMoptions.helper_funcdir,'test4switch.m');
addons{2} = fullfile(PMoptions.helper_funcdir,'timepoint_current.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{addons{1},addons{2}})
%---run-----------------------
PM_switching_model(PMoptions)
%---cleanup-------------------
delete(gcp('nocreate'))

