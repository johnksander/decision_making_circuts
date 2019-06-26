clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 3;

%---setup---------------------
tmax = 20; %diagnostics_fullnoise
options = set_options('modeltype','diagnostics','comp_location','woodstock',...
    'sim_name','test','jobID',jobID,'tmax',tmax,'cut_leave_state',tmax);

options.EtoE = .0405; %fixed 
options.ItoE = 5.8175; options.EtoI = 0.3958;
options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'
Rstim = 0; %rate for stimulus input spikes
options.trial_stimuli = [Rstim,Rstim];

%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



%when you want this code again 
% do_config = mod(options.jobID,10);
% do_config(do_config == 0) = 10;
% options.EtoE = .0405; %fixed
% options = get_network_params(do_config,options);

%options.trial_stimuli(2) = options.trial_stimuli(2) * 2;
% %adjust stimulus B strength
% stim_mod = .5:.25:2; %just randomly sample mod weight, do enough it'll even out 
% stim_mod = randsample(stim_mod,1);
% options.trial_stimuli(2) = options.trial_stimuli(2) * stim_mod; %adjust stim B
