 clear 
clc
format compact 



%---setup---------------------
jID = str2double([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 75; %trial simulation time (s) 
options = set_options('modeltype','NETS','comp_location','hpc',...
    'sim_name','nets_D2t-slower_spikedata','jobID',jID,'tmax',t,...
    'netpair_file','D2t-slower','noswitch_timeout',t,'record_spiking','on');



%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,'spikeout_model.m')
end
delete(options.output_log) %no need for these right now
% logdir = fullfile(options.save_dir,'logs'); %put them seperately
% if ~isdir(logdir),mkdir(logdir);end
%movefile(options.output_log,logdir)

% %if you need more states for specific things, use code here:
% %---need to get more states for slow nets here, new stim range---
% switch options.stim_targs
%     case 'Estay' %set to random slow network instead
%         slow_nets = 1:2:9;
%         do_config = slow_nets(randi(numel(slow_nets)));
%         options = get_network_params(do_config,options);
% end
% stim_mod = exp(1:.5:3); %new range for stim B
% stim_mod = randsample(stim_mod,1);
% options.trial_stimuli(2) = options.trial_stimuli(2) * stim_mod; %adjust stim B
% %-----------------------------

