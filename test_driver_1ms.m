clear
clc
format compact


%---setup---------------------
idx = 3; %use this to index the different network types
t = 5e3; %trial simulation time (s)
options = set_options('modeltype','equate_stim','comp_location','woodstock',...
    'timestep',2e-3,...
    'sim_name','test_real_durations','jobID',idx,'tmax',t,...
    'percent_Dslow',.5,'netpair_file','slowD',...
    'stim_pulse',[t,0],'sample_Estay_offset',0);

Rstim = 0;
options.trial_stimuli = [Rstim,Rstim];
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,'spikeout_model.m')
end
