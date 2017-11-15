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
config_options.sim_name = 'fastswitch_baseline';
options = set_options(config_options);

%get results
output_fn = [options.modeltype '_' options.sim_name '.mat'];
sim_output = load(fullfile(options.save_dir,output_fn));
options = sim_output.options; %reset to simulation options
options = reset_options_paths(options); %fix paths if needed
sim_output = sim_output.sim_results;
options.trial_hists = trial_hists;
options.conv2secs = timestep;

%find stim A & stim B durations
sim_output = cellfun(@(x) x{1},sim_output,'UniformOutput', false); %take only stay durations

%don't need to remove the first artifically induced switch anymore, because
%they're not artificially induced! 

% %get stim A durations
% stimA = cellfun(@(x) x(1:2:end),sim_output,'UniformOutput', false);
% %REMOVE THE FIRST artifically induced switch
% stimA = cellfun(@(x) x(2:end),stimA,'UniformOutput', false); %REMOVE THE FIRST artifically induced switch
% %get stim B durations
% stimB = cellfun(@(x) x(2:2:end),sim_output,'UniformOutput', false);

keyboard

%aggregate results
trial_stimvals = unique(options.trial_stimuli(:,1));
num_trials = numel(trial_stimvals);
num_stims = numel(options.trial_stimuli(1,:));
options.num_stims = num_stims;
options.stim_labels = stim_labels; %just so it's easier to pass
trial_results = cell(num_trials,num_stims);
for idx = 1:num_trials
    trials = trial_stimvals(idx) == options.trial_stimuli(:,1);
    trials = sim_output(trials);
    for stimidx = 1:numel(stim_labels)
        curr_stim = strrep(stim_labels{stimidx},'stim ','');
        trial_data = cellfun(@(x) x(strcmpi(x(:,2),curr_stim),1),trials,'UniformOutput',false);
        trial_data = cell2mat(cat(1,trial_data{:}));
        trial_results{idx,stimidx} = trial_data;
    end
    %trial_results{idx,1} = cell2mat(stimA(trials));
    %trial_results{idx,2} = cell2mat(stimB(trials));
end


%reset this shit
real_trial_currents = NaN(num_trials,num_stims);
for idx = 1:num_trials
    trials = trial_stimvals(idx) == options.trial_stimuli(:,1);
    real_trial_currents(idx,1) = unique(options.trial_stimuli(trials,1));
    real_trial_currents(idx,2) = unique(options.trial_stimuli(trials,2));
end
options.trial_currents = real_trial_currents;

%look at results
stim_duration_means = NaN(num_trials,num_stims);
stim_num_states = NaN(num_trials,num_stims);

for trialidx = 1:num_trials
    curr_trial = trial_results(trialidx,:);
    trial_histograms(options,curr_trial,trialidx);
    stim_duration_means(trialidx,:) = cellfun(@mean,curr_trial);
    stim_num_states(trialidx,:) = cellfun(@numel,curr_trial);
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


%plot it another way

close all
orient portrait 
 
fontsz = 16;
orange = [255 128 0] ./ 255;;
blue = [51 51 255] ./ 255;
plot_counter = 0;
stim_colors = {orange,blue};
for trialidx = 1:num_trials
      tile4plot = {sprintf('%s current = %.2f',options.stim_labels{1},options.trial_currents(trialidx,1)),...
         sprintf('%s current = %.2f',options.stim_labels{2},options.trial_currents(trialidx,2))}; 
    for idx = 1:num_stims
        plot_counter = plot_counter + 1;
        stim_current = options.trial_currents(trialidx,idx);
        ax(plot_counter) =  subplot(num_trials,num_stims,plot_counter);
        histogram(trial_results{trialidx,idx} * options.conv2secs,...
            'FaceColor',stim_colors{idx},'EdgeColor',stim_colors{idx})
        xlabel('state duration (s)','Fontsize',fontsz)
        ylabel('frequency','Fontsize',fontsz)
        if plot_counter == 1
            xlim([0 30]) %annoying... 
        end
        title(tile4plot,'Fontsize',fontsz)
        set(ax(plot_counter),'fontsize',fontsz)
        legend(options.stim_labels{idx})
    end
end

fig_name = 'duration_hists';
print(fullfile(options.save_dir,fig_name),'-djpeg')



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



