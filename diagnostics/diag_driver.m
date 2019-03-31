clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 1;

%---setup---------------------
tmax = 30; %diagnostics_fullnoise
options = set_options('modeltype','diagnostics','comp_location','woodstock',...
    'timestep',.02e-3,...
    'sim_name','diag_voltage_dt02','jobID',jobID,'tmax',tmax,...
    'netpair_file','fastD',...
    'stim_pulse',[tmax,0],'sample_Estay_offset',0);

do_config = mod(options.jobID,10);
do_config(do_config == 0) = 10;
options.EtoE = .0405; %fixed
options = get_network_params(do_config,options);

Rstim = 0; %rate for stimulus input spikes
options.trial_stimuli = [Rstim,Rstim];

% options.EtoE = .0405 *1;   %1
% options.ItoE = 2;%1.2904;% * 3; 
% options.EtoI = .4;%0.1948;% * 2; 
% options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'


%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



