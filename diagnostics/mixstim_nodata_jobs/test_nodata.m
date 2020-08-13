clear
clc
format compact
hold off;close all
%investigating model behavior


addpath('../')
jobID = 10;
Sname = 'test_nodata';

%---setup---------------------
tmax = 25;
options = set_options('modeltype','diagnostics','comp_location','bender',...
    'sim_name',Sname,'jobID',jobID,'tmax',tmax,'netpair_file','D2t-slower',...
    'noswitch_timeout',tmax);

%------test with stim found for network #X
do_net = 10; %fast net 5
options = get_network_params(do_net,options);
options.EtoE = .0405; %fixed

do_totals = 2; %[.5,1,2]; %do 50%, 100%, 200% total stimulus intensity 
%stim_mod = randsample(do_totals,1); %set total intensity to (stim_mod * Rstim)
stim_mod = do_totals;

%both stims are now mixed ratio 
stim_mix = {'Estay','Eswitch'}; %both targets 
p_new = .2; %proportion alternate (new) stimulus
add_targ = ~strcmp(stim_mix,options.stim_targs{1});
options.stim_targs{2} = stim_mix{add_targ};
total_strength = options.trial_stimuli{1};
total_strength = total_strength .* stim_mod;
options.trial_stimuli{1} = total_strength .* (1-p_new);
options.trial_stimuli{2} = total_strength .* p_new;


%options.trial_stimuli{1} = [0,0];
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model_lite(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now

setenv('JID',num2str(jobID))
setenv('SIM_NAME',Sname); %'diag_EtoIfixed'
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
