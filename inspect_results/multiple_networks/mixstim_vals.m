clear 
clc
format compact 



basedir = '~/Desktop/work/ACClab/rotation/project';
addpath(basedir)

%my model 
%---setup---------------------
%jID = str2double([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
jID = 2;
t = 10; %trial simulation time (s) 

options = set_options('modeltype','NETS','comp_location','bender',...
    'sim_name','test-example','jobID',jID,'tmax',t,...
    'netpair_file','D2t-slower','noswitch_timeout',t);

%%start with this network pair first...
do_net = mod(options.jobID,10);
do_net(do_net == 0) = 10;



stim_mix = {'Estay','Eswitch'}; %both targets
stims = {'A','B'};
stim_vals = array2table(NaN(0,numel(stim_mix)),'VariableNames',stims); %{'Estay','Eswitch'}


do_totals = [.5,1,2]; %do 50%, 100%, 200% total stimulus intensity
for i = 1:numel(do_totals)
    stim_mod = do_totals(i); %set total intensity to (stim_mod * Rstim)
    
    %set stimulus to mixed ratio
    mix_vals = 0:.05:1;
    for j = 1:numel(mix_vals)
        
        options = get_network_params(do_net,options);
        
        p_new = mix_vals(j); %proportion alternate (new) stimulus
        add_targ = ~strcmp(stim_mix,options.stim_targs{1});
        options.stim_targs{2} = stim_mix{add_targ};
        total_strength = options.trial_stimuli{1};
        total_strength = total_strength .* stim_mod;
        options.trial_stimuli{1} = total_strength .* (1-p_new);
        options.trial_stimuli{2} = total_strength .* p_new;
        
        curr_vals = cat(1,options.trial_stimuli{:});
        if sum(diff(curr_vals')) > 0,error('look here');end
        curr_vals = curr_vals(:,1)';
        
        curr_vals = array2table(curr_vals,'VariableNames',stims);
        stim_vals = [stim_vals;curr_vals];   
    end
    
end

AB = table2array(stim_vals);

stim_vals.total = sum(AB,2);
stim_vals.difference = abs(stim_vals.A - stim_vals.B);
stim_vals.ratio = min(AB,[],2) ./ max(AB,[],2);

fz = 16;
close all
figure;
orient tall;
subplot(2,1,1)
scatter(stim_vals.A,stim_vals.B,'filled');
ylabel('hedonic (Hz)');xlabel('aversive (Hz)')
set(gca,'FontSize',fz,'FontWeight','b')
axis square
subplot(2,1,2)
scatter(stim_vals.ratio,stim_vals.difference,'filled');
ylabel('stimulus difference');xlabel('stimulus ratio')
set(gca,'FontSize',fz,'FontWeight','b')
axis square
set(gcf,'Renderer','painters')
print('stim_info','-djpeg','-r400')


% 
% %put some info in log file 
% msg = sprintf('---stimuli mixed as %i/%i',round((1-p_new)*100),round(p_new*100));
% update_logfile(msg,options.output_log)
% for idx = 1:numel(options.stim_targs)
%     msg = sprintf('---\t %.2f Hz -> %s',options.trial_stimuli{idx}(1),options.stim_targs{idx});
%     update_logfile(msg,options.output_log)
% end
% 
% 
% %you have to checkpoint these b/c the sims are so long 
% options.grid_index = str2double(getenv('SLURM_ARRAY_TASK_ID'));%HPCC only lets indicies up to 10k!!
% checkpointFN = fullfile(options.save_dir,sprintf('checkpoint_%i.mat',options.grid_index));
% %if exist(checkpointFN) == 0 %only run unfinished jobs 
% %    delete(options.output_log);return
% %end
% 
% %---run-----------------------
% modelfile = spikeout_model_grid(options);
% %---cleanup-------------------
% if isempty(dir(fullfile(options.save_dir,'code4*zip')))
%     driverfile = mfilename;
%     backup_jobcode(options,driverfile,'spikeout_model.m')
% end
% delete(checkpointFN)
% update_logfile('checkpoint data deleted',options.output_log)
% %delete(options.output_log) %no need for these right now
% logdir = fullfile(options.save_dir,'logs'); %put them seperately
% if ~isdir(logdir),mkdir(logdir);end
% movefile(options.output_log,logdir)

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

