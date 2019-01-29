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
    'sim_name','daig_rates','jobID',jobID,'tmax',tmax,'percent_Dslow',0,...
    'stim_pulse',[tmax,0],'cut_leave_state',tmax,'sample_Estay_offset',0,...
    'ratelim_check','on');
% options = set_options('modeltype','PS','comp_location','bender',...
%     'sim_name','test_model','tmax',tmax,...
%     'stim_pulse',[tmax,0],'cut_leave_state',tmax,'sample_Estay_offset',0);


options.EtoE = .0405 *1;   %1
options.ItoE = 2.0595;%1.2904;% * 3; 
options.EtoI = .108;%0.1948;% * 2; 
options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'
Rext = 1400; %poisson spike train rate for noise, Hz
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



