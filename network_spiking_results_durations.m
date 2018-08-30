clear
clc
format compact
hold off;close all


%opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
%opt.pulse_stim = 'rem'; %'yes' | 'total_time' | 'rem' whether to treat durations as samples (rem = time during sample)

%specify simulation
%---sim setup-----------------
% Snames = {'network_spiking_P1_1','network_spiking_P2_1','network_spiking_P4_1',...
%     'network_spiking_P6_1', 'network_spiking_P8_1', 'network_spiking_P10_1',...
%     'network_spiking_P150_1'};

Snames = {'sim_v2_P1_1','sim_v2_P2_1','sim_v2_P4_1',...
    'sim_v2_P6_1','sim_v2_P8_1'};

figdir = 'network_spiking_pulse_figures';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

num_workers = numel(Snames);
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf)

parfor i = 1:numel(Snames)
    opt = struct();
    %sample duration
    opt.outcome_stat = 'mu'; opt.pulse_stim = 'yes';
    make_my_figs(basedir,Snames{i},figdir,opt)
    %duration after stim onset
    opt.outcome_stat = 'mu'; opt.pulse_stim = 'rem';
    make_my_figs(basedir,Snames{i},figdir,opt)
    %log duration after onset
    opt.outcome_stat = 'logmu'; opt.pulse_stim = 'rem';
    make_my_figs(basedir,Snames{i},figdir,opt)
    %total time
    opt.outcome_stat = 'mu'; opt.pulse_stim = 'total_time';
    make_my_figs(basedir,Snames{i},figdir,opt)
    %log total time
    opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
    make_my_figs(basedir,Snames{i},figdir,opt)
end

delete(gcp('nocreate'))

% opt = struct();
% 
% %control, log total time @ 150 pulse
% opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
% make_my_figs(basedir,'network_spiking_P150_1',figdir,opt)
% 
% %control, log total time @ 10 pulse
% opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
% make_my_figs(basedir,'network_spiking_P10_1',figdir,opt)
% 

%equate X axes for all figs of the same type
%Snames = Snames(1:end-1);
figdir =  fullfile(basedir,'Results',figdir,'durations');
savedir = fullfile(figdir,'Xsynced_dectiming');
if ~isdir(savedir),mkdir(savedir);end
ftypes = {'decision_timing_log','decision_timing','total_time','total_time_log','total_samples'};


for idx = 1:numel(ftypes)

    Xls = NaN(numel(Snames),2);
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        Xls(fidx,:) = h.CurrentAxes.XLim;
        close all
    end
    xlims = [min(Xls(:,1)),max(Xls(:,2))];
    %xlims = [0,.8];
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        ax = findobj(h,'Type','Axes');      
        set(ax,'XLim',xlims)
        
        ch = allchild(ax);
        ch = cat(1,ch{:});
        set(ch,'Normalization','probability')

        set(ch,'BinWidth',.01)
        
        print(fullfile(savedir,[Snames{fidx} '_' ftypes{idx}]),'-djpeg')
        close all
    end

end

%
% for i = 1:numel(Snames)
%     %sample duration
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'yes';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %duration after stim onset
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'rem';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %total time
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'total_time';
%     make_my_figs(basedir,Snames{i},figdir,opt)
% end


% %for summary statistics
% fprintf('network_spiking_P150_1 log(s) decision-timing\n')
% %control, log total time @ 150 pulse
% opt.outcome_stat = 'logmu'; opt.pulse_stim = 'rem'; opt.summary_stats = 'yes';
% make_my_figs(basedir,'network_spiking_P150_1',figdir,opt)
%
% fprintf('network_spiking_P101 log(s) decision-timing\n')
% %control, log total time @ 10 pulse
% opt.outcome_stat = 'logmu'; opt.pulse_stim = 'rem';  opt.summary_stats = 'yes';
% make_my_figs(basedir,'network_spiking_P10_1',figdir,opt)


function make_my_figs(basedir,sim_name,figdir,opt)
hold off;close all

outcome_stat = opt.outcome_stat;
pulse_stim = opt.pulse_stim;
print_anything = 'yes'; %'yes' | 'no';
save_anything = 'yes';
summary_stats = 'no'; %summary_stats = opt.summary_stats;

figdir = fullfile(basedir,'Results',figdir,'durations');
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_BL'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);

%pulse duration... kinda hardcoded here
Pdur = strsplit(sim_name,'_');
Pdur = Pdur{3};
Pdur = strrep(Pdur,'P','');
Pdur = sprintf('%ss',Pdur);


timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
%this file was being created "by hand" with parsweep_find_examples
%network_pair_info = load(fullfile(basedir,'helper_functions','network_pairs'));
%network_pair_info = network_pair_info.network_pairs;
%take this stuff directly from get_network_params()
%look in get_network_params() for this NPjobs, 0 is the last one paired w/ 9..
NPjobs = cat(1,[1:2:9],[2:2:9,0])'; %u-g-l-y
network_pair_info = cell(size(NPjobs,1),1);
for idx = 1:numel(NPjobs)/2
    CP = num2cell(NPjobs(idx,:))';
    CP = cellfun(@(x,y) {get_network_params(x,y)}, CP,repmat({struct()},2,1));
    CP = cellfun(@(x) {x.ItoE,x.EtoI,unique(x.trial_stimuli),x.stim_targs},...
        CP,'UniformOutput',false);
    network_pair_info{idx} = cat(1,CP{:});
end

%this code is holdover from when it was loaded.. also dumb.
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],...
    network_pair_info,'UniformOutput',false);
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],...
    network_pair_info,'UniformOutput',false);

fprintf('\n---loading simulation: %s\n',sim_name)

%get results
file_data = cell(num_files,2);
for idx = 1:num_files
    if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %store parameters
    file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{1};
        
    %     %--this is prior to change in durations data structure
    %     %take only stimulus state durations
    %     state_durations = state_durations{1}(:,1);
    %     %     state_durations = cellfun(@(x) x(:,1),state_durations,'UniformOutput',false);
    %     %     state_durations = vertcat(state_durations{:});
    %     state_durations = cat(1,state_durations{:}); %ooo that's annoying
    %
    %     state_durations = find_stay_durations(state_durations,curr_file.options);
    %
    %     %convert to time
    %     state_durations = state_durations * timestep;
    %     switch pulse_stim
    %         case 'yes' %just do this now while options is handy
    %             state_durations = floor(state_durations ./ sum(curr_file.options.stim_pulse));
    %         case 'rem' %look at when IN the sample switch happened
    %             state_durations = mod(state_durations,sum(curr_file.options.stim_pulse));
    %         case 'total_time' %need to subtract the off-period time (in betwn pulses)
    %             Nsamps = floor(state_durations ./ sum(curr_file.options.stim_pulse));
    %             state_durations = state_durations - (Nsamps.*curr_file.options.stim_pulse(2));
    %     end

    state_durations = find_stay_durations(state_durations,curr_file.options);
    switch pulse_stim
        case 'yes' %just do this now while options is handy
            state_durations = state_durations.samples;
        case 'rem' %look at when IN the sample switch happened
            state_durations = state_durations.decision_time;
        case 'total_time'
            state_durations = state_durations.duration;
    end
    
    file_data{idx,1} = state_durations;
end

%search for jobs with identical parameters, collapse distributions
%get the randomized network parameters
job_params = cellfun(@(x)...
    [x.ItoE, x.EtoI,unique(x.trial_stimuli),find(strcmpi(x.stim_targs, stimtarg_vals))],...
    file_data(:,2),'UniformOutput',false); %matching "network_pair_info" format
job_params = vertcat(job_params{:});
uniq_params = unique(job_params,'rows');
num_jobs = size(uniq_params,1);
fprintf('----------------------\n')
fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)

%collapse duplicate job parameters
result_data = cell(num_jobs,2);
Nruns = NaN(num_jobs,1); %record the number of successful jobs..
for idx = 1:num_jobs
    %find all matching
    curr_file = ismember(job_params,uniq_params(idx,:),'rows');
    Nruns(idx) = sum(curr_file);
    fprintf('---\nparameter set %.3f %.3f %.3f %s, n files = %i\n',uniq_params(idx,1:3),...
        stimtarg_vals{uniq_params(idx,4)},Nruns(idx))
    %collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    fprintf('------n states = %i\n',numel(result_data{idx,1}))
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
end

% %just grab some simple stats here
% mu_duration = cellfun(@(x) mean(x),result_data(:,1));
% med_duration = cellfun(@(x) median(x),result_data(:,1));
% logmu_dur = cellfun(@(x) mean(log(x)),result_data(:,1));
% %get the network parameters
% EtoE = cellfun(@(x)  x.EtoE,result_data(:,2));
% ItoE = cellfun(@(x)  x.ItoE,result_data(:,2));
% EtoI = cellfun(@(x)  x.EtoI,result_data(:,2));
% stim_value = cellfun(@(x)  unique(x.trial_stimuli),result_data(:,2));

net_type = num2cell(num2cell(uniq_params),2);
net_type = cellfun(@(x)  [x(1:3),stimtarg_vals{x{4}}],net_type,'UniformOutput',false);
net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],net_type,'UniformOutput',false);
net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],net_type,'UniformOutput',false);

% fig_fn = 'duration_dists';
%
% switch pulse_stim
%     case 'rem' %rename for special case
%         fig_fn = 'duration_after_stim_onset';
% end

fig_fn = [sim_name '_%s'];

switch pulse_stim
    case 'yes'
        unit_measure = 'samples';
        fig_fn = sprintf(fig_fn,'total_samples');
    case 'rem'
        unit_measure = 's - onset';
        fig_fn = sprintf(fig_fn,'decision_timing');
    otherwise
        unit_measure = 's'; %like "seconds" not samples
        fig_fn = sprintf(fig_fn,'total_time');
end


switch outcome_stat
    case 'mu'
        %outcome = mu_duration;
        %Zlabel = 'mean duration (s)';
        Zlabel = sprintf('state duration (%s)',unit_measure);
        %figdir = fullfile(figdir,'mean_duration');
    case 'med'
        %outcome = med_duration;
        Zlabel =  sprintf('median duration (%s)',unit_measure);
        %figdir = fullfile(figdir,'med_duration');
    case 'logmu'
        %outcome = logmu_dur;
        %Zlabel = 'mean log duration [ log(s) ]';
        Zlabel = sprintf('log(%s) state duration',unit_measure);
        %figdir = fullfile(figdir,'logmean_duration');
        fig_fn = [fig_fn,'_log'];
end

if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;
alph = .5;

num_types = numel(network_pair_info);
plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        hold on
        
        %get the right color
        if strcmpi(curr_net_info{j,end},'slow')
            lcol = matorange;
        elseif strcmpi(curr_net_info{j,end},'fast')
            lcol = matblue;
        end
        
        %find the right results for network set-up
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = cell2mat(result_data(curr_data,1));
        
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = curr_data(curr_data ~= 0);%inf errors
                    curr_data = log(curr_data);
            end
            histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
            Ylab = 'freq';
            
            %             switch pulse_stim
            %                 case {'yes','rem'}
            %                     histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
            %                     Ylab = 'freq';
            %                 otherwise
            %                     [kde,kde_i] = ksdensity(curr_data);
            %                     area(kde_i,kde,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
            %                     Ylab = 'p(x)';
            %             end
        end
        %take control data
        BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
        curr_data = cellfun(@(x) isequal(x,BLinfo),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = cell2mat(result_data(curr_data,1));
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = curr_data(curr_data ~= 0);%inf errors
                    curr_data = log(curr_data);
            end
            switch pulse_stim
                case 'yes'
                    histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
                    Ylab = 'freq';
                otherwise
                    [kde,kde_i] = ksdensity(curr_data);
                    area(kde_i,kde,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
                    Ylab = 'p(x)';
            end
        end
        hold off
        legend(sprintf('%.0fHz %s',curr_net_info{j,3},Pdur),'location','best')
        
        if plt_idx == 9 || plt_idx == 10
            xlabel(Zlabel)
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info{j,4}),'Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('network #%i %s',idx,Ylab))
        end
        
    end
end
orient tall
linkaxes(h,'x')
axis tight

switch print_anything
    case 'yes'
        print(fullfile(figdir,fig_fn),'-djpeg')
        savefig(fullfile(figdir,fig_fn))
end

switch save_anything
    case 'yes'
        svdir = fullfile(figdir,'data');
        if ~isdir(svdir),mkdir(svdir);end
        svFN = [sim_name '_%s'];
        switch pulse_stim
            case 'yes'
                svFN = sprintf(svFN,'total_samples');
            case 'rem'
                svFN = sprintf(svFN,'decision_timing');
            otherwise
                svFN = sprintf(svFN,'total_time');
        end
        save(fullfile(svdir,svFN),'result_data','network_pair_info')
        
end
%info about simulation
num_states = cellfun(@(x) numel(x{1}),num2cell(result_data,2));
need_more = num_states < 10000;
need_more = net_type(need_more);
current_count = num_states(num_states < 10000);
for idx = 1:num_types
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),need_more,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        if sum(curr_data) > 0
            fprintf('\nnetwork #%i %s has < 10k states (%i)',idx,curr_net_info{j,4},current_count(curr_data))
            fprintf('\n---parameter set:\n')
            disp(curr_net_info(j,:))
            fprintf('\n------------------------\n')
        end
        
        BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
        curr_data = cellfun(@(x) isequal(x,BLinfo),need_more,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        if sum(curr_data) > 0
            fprintf('\nnetwork #%i %s has < 10k states (%i)',idx,BLinfo{4},current_count(curr_data))
            fprintf('\n---parameter set:\n')
            disp(BLinfo)
            fprintf('\n------------------------\n')
        end
        
        
    end
    
end

switch summary_stats
    case 'yes'
        %summary stats
        fprintf('\n------------------------\n')
        fprintf('Summary statistics\n')
        fprintf('%s:\n',Zlabel)
        fprintf('------------------------\n\n')
        for idx = 1:num_types
            fprintf('\n------------------------\nNetwork #%i\n',idx)
            curr_net_info = network_pair_info{idx};
            Xstim = cell(2,1); %for testing difference between stim distributions
            for j = 1:2
                fprintf('---type: %s\n',curr_net_info{j,end})
                %find the right results for network set-up
                curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
                curr_data = cat(1,curr_data{:});
                curr_data = result_data{curr_data,1};
                switch outcome_stat
                    case 'logmu'
                        curr_data = log(curr_data(curr_data~=0));
                end
                Xstim{j} = curr_data;
                
                fprintf('            stimulus = %.2f\n',curr_net_info{j,3})
                fprintf('            ---mean = %.2f\n',mean(curr_data))
                
                
                %                 %take control data
                %                 BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
                %                 curr_data = cellfun(@(x) isequal(x,BLinfo),net_type,'UniformOutput',false);
                %                 curr_data = cat(1,curr_data{:});
                %                 fprintf('            baseline = %.2f\n',outcome(curr_data))
            end
            
            fprintf('\n---hyp. test: mu stim durrations\n')
            switch outcome_stat
                case 'logmu'
                    Xstim = cellfun(@(x) log(x(x~=0)),Xstim,'UniformOutput',false);
            end
            [~,pval] = ttest2(Xstim{1},Xstim{2}); %regular old t-test
            fprintf('t-test p = %.3f\n',pval)
            [CI,H] = boot_mean_diffs(Xstim{1},Xstim{2},10000);
            fprintf('bootstrap test: %s\n',H)
            fprintf('bootstrap CI: %.2f, %.2f\n',CI)
        end
end
end


