clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 1;

%---setup---------------------
tmax = 20;
options = set_options('modeltype','diagnostics','comp_location','bender',...
    'sim_name','pulse_test','jobID',jobID,'tmax',tmax,'netpair_file','D2t-slower',...
    'stim_pulse',[1,3]);

%------test with stim found for network #1 
options = get_network_params(1,options);
options.EtoE = .0405; %fixed
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now

setenv('JID',num2str(jobID))
setenv('SIM_NAME',options.sim_name); %'diag_EtoIfixed'
inspect

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
