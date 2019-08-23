clear 
clc
format compact 




%my model
%---setup---------------------
jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 350; %trial simulation time (s)
options = set_options('modeltype','NETS','comp_location','hpc',...
    'sim_name','nets_fastD_baseline','jobID',jID,'tmax',t,...
    'netpair_file','fastD',...
    'record_spiking','on');

%specify baseline stimulus 0hz
options.stim_targs = 'baseline'; 
Rstim = 0; 
options.trial_stimuli = [Rstim,Rstim];

%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
delete(options.output_log) %no need for these right now
% logdir = fullfile(options.save_dir,'logs'); %put them seperately
% if ~isdir(logdir),mkdir(logdir);end
% movefile(options.output_log,logdir)