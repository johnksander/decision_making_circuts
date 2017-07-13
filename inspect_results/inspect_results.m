clear
clc
format compact

%result summaries
fontsz = 16;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%specify simulation
%---sim setup-----------------
config_options.modeltype = 'JK';
config_options.sim_name = 'switching_bias8_reparam2';
options = set_options(config_options);

%get results
output_fn = [options.modeltype '_' options.sim_name '.mat'];
sim_output = load(fullfile(options.save_dir,output_fn));
options = sim_output.options; %reset to simulation options
options = reset_options_paths(options); %fix paths if needed
sim_output = sim_output.sim_results;
options.trial_hists = trial_hists;
options.conv2secs = timestep;

%look at results
num_trials = numel(options.trial_currents(:,1));
num_stims = numel(options.trial_currents(1,:));
options.num_stims = num_stims;
options.stim_labels = stim_labels; %just so it's easier to pass
stim_duration_means = NaN(num_trials,num_stims);
stim_num_states = NaN(num_trials,num_stims);

for trialidx = 1:num_trials
    
    trial_results = sim_output{trialidx};
    trial_results = trial_results{1}; %just take the "stay" state
    
    %trial_results = flipud(trial_results); %you're an idiot.
    
    %sort by stimuli
    Ainds = 1:2:numel(trial_results);
    Ainds = Ainds(2:end); %REMOVE the first artifically induced switch 
    Binds = 2:2:numel(trial_results);
    trial_results = {trial_results(Ainds),trial_results(Binds)};
    %do functions
    trial_histograms(options,trial_results,trialidx);
    stim_duration_means(trialidx,:) = cellfun(@mean,trial_results);
    stim_num_states(trialidx,:) = cellfun(@numel,trial_results);
end
close all

stim_duration_means = stim_duration_means * options.conv2secs; %take care of unit conversion here
stim_duration_means(isnan(stim_duration_means)) = 0;
stim_num_states(isnan(stim_num_states)) = 0;
num_switches = sum(stim_num_states,2) - 1;
num_switches(num_switches == -1) = 0; %if there were no switches...

figure(1)

curr_diff = options.trial_currents(:,1) - options.trial_currents(:,2);
for stimidx = 1:num_stims
    plot(curr_diff,stim_duration_means(:,stimidx),'linewidth',2)
    hold on 
end
% for stimidx = 1:num_stims
%     plot(options.trial_currents(:,stimidx),stim_duration_means(:,stimidx),'linewidth',2)
%     hold on %!MUST FIX FOR REAL SIMULATION ^^^
% end
ylabel('Average state duration (s)')
xlabel('Stim A current bias')
title('State durations')
legend(stim_labels)
set(gca,'Fontsize',fontsz)
orient landscape
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.0f%%',x*100),Xticks,'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
hold off
fig_fn = [options.sim_name '_mean_durations'];
print(fullfile(options.save_dir,fig_fn),'-djpeg')



figure(2)
plot(curr_diff,num_switches,'linewidth',2) 
title('Stimuli switching')
ylabel('number of switches between stimuli')
xlabel('Stim A current bias')
set(gca,'Fontsize',fontsz)
orient landscape
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.0f%%',x*100),Xticks,'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
fig_fn = [options.sim_name '_num_switches'];
print(fullfile(options.save_dir,fig_fn),'-djpeg')


% figure(2)
% plot(options.trial_currents(:,1),num_switches,'linewidth',2) %!MUST FIX FOR REAL SIMULATION
% title('State switching')
% ylabel('number of state switches')
% xlabel('Current (proportion)')
% set(gca,'Fontsize',fontsz)
% orient landscape
% fig_fn = [options.sim_name '_num_switches'];
% print(fullfile(options.save_dir,fig_fn),'-djpeg')


% figure(3)
% 
% curr_diff = options.trial_currents(:,1) - options.trial_currents(:,2);
% for stimidx = 1:num_stims
%     plot(curr_diff,stim_duration_means(:,stimidx),'linewidth',2)
%     hold on 
% end
% ylim([0 100])
% ylabel('Average state duration (s)')
% xlabel('Stim A current bias')
% title('State durations')
% legend(stim_labels)
% set(gca,'Fontsize',fontsz)
% orient landscape
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%.0f%%',x*100),Xticks,'UniformOutput', false);
% set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
% hold off
% fig_fn = [options.sim_name '_mean_durationsYlim'];
% print(fullfile(options.save_dir,fig_fn),'-djpeg')



