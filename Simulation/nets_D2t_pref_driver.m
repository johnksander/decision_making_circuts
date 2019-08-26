clear 
clc
format compact 




%my model 
%---setup---------------------
jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 350; %trial simulation time (s) 
options = set_options('modeltype','NETS','comp_location','hpc',...
    'sim_name','nets_D2t_20s_pref','jobID',jID,'tmax',t,...
    'netpair_file','D2t_20s','record_spiking','off');


%adjust stimulus B strength
stim_mod = 0:.25:2; % 0:.25:2; %just randomly sample mod weight, do enough it'll even out 
stim_mod = randsample(stim_mod,1);
options.trial_stimuli(2) = options.trial_stimuli(2) * stim_mod; %adjust stim B
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
delete(options.output_log) %no need for these right now
