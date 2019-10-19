clear 
clc
format compact 

profile on

%my model 
%---setup---------------------
% %jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
% t = 30; %trial simulation time (s) 
% options = set_options('modeltype','PS','comp_location','woodstock',...
%     'sim_name','profile_test','tmax',t,...
%     'ratelim_check','on','cut_leave_state',t);

%my model 
%---setup---------------------
t = 10; %trial simulation time (s) 
options = set_options('modeltype','NETS','comp_location','woodstock',...
    'sim_name','profile_test','tmax',t,...
    'netpair_file','D2t','record_spiking','off');

%------test with stim found for network #1 
options = get_network_params(1,options);
Rstim = 68; %found in search 
options.trial_stimuli = [Rstim,Rstim];


%---run-----------------------
%modelfile = spikeout_model(options);
%modelfile = spikeout_model_GPU(options);
modelfile = spikeout_model_GPU_single(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)


profile off
profsave(profile('info'),'MYPROFILE_RESULTS')
