clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 1;

%---setup---------------------
tmax = 60; %diagnostics_fullnoise
options = set_options('modeltype','diagnostics','comp_location','woodstock',...
    'sim_name','diag_pulses','jobID',jobID,'tmax',tmax,...
    'percent_Dslow',.5,'netpair_file','slowD',...
    'state_def','include_undecided','stim_pulse',[2,1]);

do_config = mod(options.jobID,10);
do_config(do_config == 0) = 10;
options.EtoE = .0405; %fixed
options = get_network_params(do_config,options);

%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



