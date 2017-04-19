clear 
clc
format compact 




%Paul's model 
%---setup---------------------
config_options.modeltype = 'PM';
config_options.sim_name = 'test2';
options = set_options(config_options);
options.equal_pools = 'off'; %'on' | 'off' set switch & stay to have equal properties (b)
options.current_pulse = 'off';  %'on' | 'off' set off for low memory implementation
options.tmax = 250 * 1000; %trial simulation time (ms) 
stimA = [1:.02:1.06]'; stimB = ones(size(stimA)); %current modifiers for stimuli (col vectors!)
options.trial_currents = [stimA,stimB]; %must be trials X stims (num_trials = numel(rows);)
options.parforlog = 'on'; 
%---parfor--------------------
num_workers = 4; %parpool workers
addons{1} = fullfile(options.helper_funcdir,'test4switch.m');
addons{2} = fullfile(options.helper_funcdir,'timepoint_current.m');
addons{3} = fullfile(options.helper_funcdir,'next_timepoint.m');
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',addons)
%---run-----------------------
modelfile = PM_model(options);
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(gcp('nocreate'))


