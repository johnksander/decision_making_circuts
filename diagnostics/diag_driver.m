clear
clc
format compact

%investigating model behavior

addpath('../')

%my model
%---setup---------------------
config_options.modeltype = 'PS_stim';
config_options.sim_name = 'diagnostics';
%specify network #1 slow w/ jobID
config_options.jobID = 2; %str2num(getenv('SGE_TASK_ID'));
config_options.tmax = 30; %trial simulation time (s)
config_options.force_back2stay = true;
config_options.stim_pulse = [1, 1]; %on, off (s)
options = set_options(config_options);

%NOTE: I AM NOT removing the first artificial stay state here. 
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now

% 
% %just look at some stuff here 
% cd(options.save_dir)
% load('PS_diagnostics_1.mat');load('PS_diagnostics_1_D.mat');
% load('PS_diagnostics_1_V.mat');load('PS_diagnostics_1_spikes.mat');
% last_idx = isnan(Vrec(1,:));
% last_idx = find(last_idx,1,'first');
% 
% timestep = .25e-3; %.25 milisecond timestep
% spikeplot = make_spikeplot(spikes(:,1:last_idx));
% imagesc(spikeplot)
% 
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% title({'spikes','(spikes in matrix enlarged for visualization)'})
% ylabel('cell')
% xlabel('time (s)')
% 
% 
% figure;
% plot(Drec(1,1:last_idx),'LineWidth',2)
% Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% xlabel('time (s)')
% 
% durations = sim_results{1};
% timecourse = cell(5,2);
% timecourse(1:2:end,:) = durations{1};
% timecourse(2:2:end,:) = durations{2};
% event_times = cell2mat(timecourse(:,1));
% event_times = cumsum(event_times) * timestep;
% timecourse = [num2cell(event_times),timecourse];






