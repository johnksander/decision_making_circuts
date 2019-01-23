clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
jobID = 3;

%---setup---------------------
tmax = 30; %diagnostics_fullnoise
options = set_options('modeltype','diagnostics','comp_location','woodstock',...
    'sim_name','diag_EtoIfixed','jobID',jobID,'tmax',tmax,...
    'stim_pulse',[tmax,0],'cut_leave_state',tmax,'sample_Estay_offset',0);

options.EtoE = .0405 *1;   %1
options.ItoE = 5; %0.885;
options.EtoI = 0.1948; %1.492
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



