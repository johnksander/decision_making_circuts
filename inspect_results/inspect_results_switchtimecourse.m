clear
clc
close all
format compact

%result summaries
fontsz = 12;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];

%specify simulation
%---sim setup-----------------
config_options.modeltype = 'JK';
config_options.sim_name = 'switching_dynamics_EIunbal';
options = set_options(config_options);

%get results
output_fn = [options.modeltype '_' options.sim_name '.mat'];
sim_output = load(fullfile(options.save_dir,output_fn));
options = sim_output.options; %reset to simulation options
options = reset_options_paths(options); %fix paths if needed
sim_timecourse_output = sim_output.sim_switch_timecourses; %get switching dynamics data
sim_output = sim_output.sim_results;
options.trial_hists = trial_hists;
options.conv2secs = timestep;

%concatenate results across runs 
summed_timecourses = cellfun(@(x) x{1}, sim_timecourse_output, 'UniformOutput', false);
summed_timecourses = cat(3,summed_timecourses{:});
summed_timecourses = sum(summed_timecourses,3);

stateswich_counts = cellfun(@(x) x{2}, sim_timecourse_output, 'UniformOutput', false);
stateswich_counts = sum(cell2mat(stateswich_counts));

avg_timecourse = summed_timecourses ./ stateswich_counts;

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 150;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [2/3 1/3]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);

%break it down by celltype & aggregate 
Estay = celltype.excit & celltype.pool_stay;
Eswitch = celltype.excit & celltype.pool_switch;
Istay = celltype.inhib & celltype.pool_stay;
Iswitch = celltype.inhib & celltype.pool_switch;

Estay = mean(avg_timecourse(Estay,:));
Eswitch = mean(avg_timecourse(Eswitch,:));
Istay = mean(avg_timecourse(Istay,:));
Iswitch = mean(avg_timecourse(Iswitch,:));

%plot the aggregated timecourses
figure(1) 
hold on
plot(Estay,'Linewidth',2)
plot(Eswitch,'Linewidth',2)
plot(Istay,'Linewidth',2)
plot(Iswitch,'Linewidth',2)
hold off
legend({'Excit. stay','Excit. switch','Inhib. stay','Inhib. switch'},'location','southwest',...
    'Fontsize',fontsz)

onset_switch = 250e-3/timestep; %add the switching time 
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
Yticks = num2cell(get(gca,'Ytick'));
Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
xlabel('state switch timecourse (ms)')
ylabel('current noise (picoamps)')
title({'cell-type mean timecourse';'switching from stimuli A (1.5 x current) to B (1 x current)'})
y_range = get(gca,'YLim');
x_range = get(gca,'XLim');
text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
set(gca,'Fontsize',fontsz)
fig_fn = 'switching_timecourses';
print(fullfile(options.save_dir,fig_fn),'-djpeg')

%now with a smoothed line
%plot the aggregated timecourses
figure(2) 
hold on
smoothfac = 10;
plot(smooth(Estay,smoothfac),'Linewidth',2)
hold on
plot(smooth(Eswitch,smoothfac),'Linewidth',2)
plot(smooth(Istay,smoothfac),'Linewidth',2)
plot(smooth(Iswitch,smoothfac),'Linewidth',2)
hold off
legend({'Excit. stay','Excit. switch','Inhib. stay','Inhib. switch'},'location','southwest',...
    'Fontsize',fontsz)

onset_switch = 250e-3/timestep; %add the switching time 
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
Yticks = num2cell(get(gca,'Ytick'));
Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
xlabel('state switch timecourse (ms)')
ylabel('current noise (picoamps)')
title({'cell-type mean timecourse';'switching from stimuli A (1.5 x current) to B (1 x current)';...
    sprintf('plot smoothed: %i-sample boxcar conv',smoothfac)})
y_range = get(gca,'YLim');
x_range = get(gca,'XLim');
text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
set(gca,'Fontsize',fontsz)
fig_fn = 'switching_timecourses_smoothed';
print(fullfile(options.save_dir,fig_fn),'-djpeg')



