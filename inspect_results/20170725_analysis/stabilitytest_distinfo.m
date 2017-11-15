clear
clc
format compact

%result summaries
fontsz = 16;
timestep = .25e-3;
home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/';
addpath(home_dir)
addpath(fullfile(home_dir,'helper_functions'))
fig_dir = fullfile(home_dir,'Results/change_switchspeed');
% stat_file = fullfile(fig_dir,'switch_info.txt');
% if ~isdir(fig_dir),mkdir(fig_dir);end
% fig_name = 'switchtime_hists';

%this is coded a little specific to the recent simulations...

%specify simulation
%---sim setup-----------------
sim_names = {'fastswitch_baseline','slowswitch_baseline',...
    'fastswitch_stimulus','slowswitch_stimulus'};
num_sims = numel(sim_names);
sim_data = load_simdata(sim_names);

orient landscape
for idx = 1:num_sims
    options = sim_data{idx,1};
    trial_data = sim_data{idx,2};
    stim_value = unique(options.trial_stimuli);
    ax(idx) =  subplot(num_sims/2,2,idx);
    histogram(trial_data * timestep)
    if mod(idx,2) == 1
        ylabel('frequency','Fontsize',fontsz)
    end
    if idx > num_sims - 2
        xlabel('state duration (s)','Fontsize',fontsz)
    end
    sim_label = sim_names{idx};
    sim_label = strrep(sim_label,'_',' ');
    tile4plot = {sim_label, [' (' num2str(stim_value)...
        ' Hz stimulation)']};
    title(tile4plot,'Fontsize',fontsz)
    set(ax(idx),'fontsize',fontsz)
    keyboard
    %print info to console
%     txtappend(stat_file,sprintf('\n------simulation: %s\n',sim_label))
%     txtappend(stat_file,'state durations\n')
%     txtappend(stat_file,sprintf('   mean = %.2f s\n',mean(trial_data * timestep)));
%     txtappend(stat_file,sprintf('   median = %.2f s\n',median(trial_data * timestep)));
%     txtappend(stat_file,'switches\n')
%     txtappend(stat_file,sprintf('   n = %i\n\n',numel(trial_data)))
end

pause(.1)




function sim_data = load_simdata(sim_names)
num_sims = numel(sim_names);
sim_data = cell(num_sims,1);

for loadidx = 1:num_sims
    
    config_options.modeltype = 'JK';
    config_options.sim_name = sim_names{loadidx};
    options = set_options(config_options);
    %get results
    output_fn = [options.modeltype '_' options.sim_name '.mat'];
    sim_output = load(fullfile(options.save_dir,output_fn));
    options = sim_output.options; %reset to simulation options
    options = reset_options_paths(options); %fix paths if needed
    sim_output = sim_output.sim_results;
    %find stim A & stim B durations
    sim_output = cellfun(@(x) x{1},sim_output,'UniformOutput', false); %take only stay durations
    sim_output = cat(1,sim_output{:});
    sim_output = cell2mat(cat(1,sim_output(:,1)));
    sim_data{loadidx,1} = options;
    sim_data{loadidx,2} = sim_output;
end

end









