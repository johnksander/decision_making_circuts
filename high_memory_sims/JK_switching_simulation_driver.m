clear 
clc
format compact 




%my model 
%---setup---------------------
config_options.modeltype = 'JK';
config_options.sim_name = 'find_baseline3';
JKoptions = set_options(config_options);
JKoptions.current_pulse = 'off';  %'on' | 'off' set off for low memory implementation
JKoptions.tmax = 1000; %trial simulation time (s) 
stimA = [.5:.02:1.5]'; stimB = [.5:.02:1.5]'; %current modifiers for stimuli (col vectors!)
JKoptions.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
JKoptions.parforlog = 'on'; 
%---parfor--------------------
num_workers = 12; %parpool workers
addons{1} = fullfile(JKoptions.helper_funcdir,'test4switch.m');
addons{2} = fullfile(JKoptions.helper_funcdir,'timepoint_current.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{addons{1},addons{2}})
%---run-----------------------
JK_switching_model(JKoptions)
%---cleanup-------------------
delete(gcp('nocreate'))

