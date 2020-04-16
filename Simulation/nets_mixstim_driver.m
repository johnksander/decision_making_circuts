clear 
clc
format compact 




%my model 
%---setup---------------------
jID = str2double([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 1500; %trial simulation time (s) 

options = set_options('modeltype','NETS','comp_location','hpc',...
    'sim_name','nets_mixstim','jobID',jID,'tmax',t,...
    'netpair_file','D2t-slower','noswitch_timeout',t);

%%start with this network pair first...
do_nets = [3,4];
options = get_network_params(randsample(do_nets,1),options);
stim_mod = .5; %set total intensity to (stim_mod * Rstim)

%set stimulus to mixed ratio 
mix_vals = 0:.05:1;
stim_mix = {'Estay','Eswitch'}; %both targets 
p_new = randsample(mix_vals,1); %proportion alternate (new) stimulus
add_targ = ~strcmp(stim_mix,options.stim_targs{1});
options.stim_targs{2} = stim_mix{add_targ};
total_strength = options.trial_stimuli{1};
total_strength = total_strength .* stim_mod;
options.trial_stimuli{1} = total_strength .* (1-p_new);
options.trial_stimuli{2} = total_strength .* p_new;

%put some info in log file 
msg = sprintf('---stimuli mixed as %i/%i',round((1-p_new)*100),round(p_new*100));
update_logfile(msg,options.output_log)
for idx = 1:numel(options.stim_targs)
    msg = sprintf('---\t %.2f Hz -> %s',options.trial_stimuli{idx}(1),options.stim_targs{idx});
    update_logfile(msg,options.output_log)
end


%you have to checkpoint these b/c the sims are so long 
options.grid_index = str2double(getenv('SLURM_ARRAY_TASK_ID'));%HPCC only lets indicies up to 10k!!
checkpointFN = fullfile(options.save_dir,sprintf('checkpoint_%i.mat',options.grid_index));
%if exist(checkpointFN) == 0 %only run unfinished jobs 
%    delete(options.output_log);return
%end

%---run-----------------------
modelfile = spikeout_model_grid(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,'spikeout_model.m')
end
delete(checkpointFN)
update_logfile('checkpoint data deleted',options.output_log)
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)

% %if you need more states for specific things, use code here:
% %---need to get more states for slow nets here, new stim range---
% switch options.stim_targs{1}
%     case 'Estay' %set to random slow network instead
%         slow_nets = 1:2:9;
%         do_config = slow_nets(randi(numel(slow_nets)));
%         options = get_network_params(do_config,options);
% end
% stim_mod = exp(1:.5:3); %new range for stim B
% stim_mod = randsample(stim_mod,1);
% options.trial_stimuli{1}(2) = options.trial_stimuli{1}(2) * stim_mod; %adjust stim B
% %-----------------------------

