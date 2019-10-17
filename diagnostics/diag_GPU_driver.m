clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 5;

%---setup---------------------
tmax = 10; %diagnostics_fullnoise
options = set_options('modeltype','diagnostics','comp_location','woodstock',...
    'sim_name','test','jobID',jobID,'tmax',tmax, 'netpair_file','D2t');

%------test with stim found for network #1 
options = get_network_params(1,options);
options.EtoE = .0405; %fixed
Rstim = 68; %found in search 
options.trial_stimuli = [Rstim,Rstim];

%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model_GPU(options);
    %[modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now

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
