clear
clc
format compact
hold off;close all


opt = struct();
opt.print_anything = 'no'; %'yes' | 'no';
opt.save_anything = 'no';
opt.multiple_stimuli = 'yes'; 
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)


%specify simulation
%---sim setup-----------------
% Snames = {'network_spiking_P1_1','network_spiking_P2_1','network_spiking_P4_1',...
%     'network_spiking_P6_1', 'network_spiking_P8_1', 'network_spiking_P10_1',...
%     'network_spiking_P150_1'};

Snames = {'dttest_0-1_fastD'};
figdir = 'figures_dttest'; 
basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

%control, log total time @ 150 pulse
opt.outcome_stat = 'mu'; opt.pulse_stim = 'off';
make_my_figs(basedir,Snames{1},figdir,opt)

% num_workers = numel(Snames);
% c = parcluster('local');
% c.NumWorkers = num_workers;
% parpool(c,c.NumWorkers,'IdleTimeout',Inf)
% 
% parfor i = 1:numel(Snames)
%     opt = struct();
%     %sample duration
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'yes';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %duration after stim onset
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'rem';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %log duration after onset
%     opt.outcome_stat = 'logmu'; opt.pulse_stim = 'rem';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %total time
%     opt.outcome_stat = 'mu'; opt.pulse_stim = 'total_time';
%     make_my_figs(basedir,Snames{i},figdir,opt)
%     %log total time
%     opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
%     make_my_figs(basedir,Snames{i},figdir,opt)
% end
% 
% delete(gcp('nocreate'))

%control, log total time @ 150 pulse
opt.outcome_stat = 'mu'; opt.pulse_stim = 'off';
make_my_figs(basedir,Snames{1},figdir,opt)

%control, log total time @ 10 pulse
opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
make_my_figs(basedir,Snames{1},figdir,opt)


%equate X axes for all figs of the same type
%Snames = Snames(1:end-1);
figdir =  fullfile(basedir,'Results',figdir,'durations');
savedir = fullfile(figdir,'Xsynced');
if ~isdir(savedir),mkdir(savedir);end
ftypes = {'decision_timing_log','decision_timing','total_time','total_time_log','total_samples'};

keyboard
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
        Fdir = fullfile(savedir,ftypes{idx});
        if ~isdir(Fdir),mkdir(Fdir);end
        print(fullfile(Fdir,[Snames{fidx} '_' ftypes{idx}]),'-djpeg')
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
print_anything = opt.print_anything; %'yes' | 'no';
save_anything = opt.save_anything;
summary_stats = 'no'; %summary_stats = opt.summary_stats;

helpdir = fullfile(basedir,'helper_functions');
figdir = fullfile(basedir,'Results',figdir,'durations');
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_BL'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);
param_varnams = {'ItoE','EtoI','stim_A','stim_B','stim_targs'};
%get general options file from the first file 
gen_options = load(output_fns{1});
gen_options = gen_options.options;
timestep = gen_options.timestep;


switch pulse_stim
    case 'off' 
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        Pdur = strsplit(sim_name,'_');
        Pdur = Pdur{3};
        Pdur = strrep(Pdur,'P','');
        Pdur = sprintf('%ss',Pdur);
end


%timestep is in set_options() now

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'}; 
keyboard
network_pair_info = load(fullfile(helpdir,'network_pairs',gen_options.netpair_file));
network_pair_info = network_pair_info.network_pairs;
num_net_types = cellfun(@(x) size(x,1),network_pair_info);
num_net_types = sum(num_net_types);
num_pairs = numel(network_pair_info);

%add a little more info to the table 
for idx = 1:num_pairs
    T = network_pair_info{idx};
    T_stargs = T.Properties.RowNames;
    T_stargs = strrep(T_stargs,'fast','Estay');
    T_stargs = strrep(T_stargs,'slow','Eswitch');
    T.stim_targs = T_stargs;
    network_pair_info{idx} = T;
end


% %this file was being created "by hand" with parsweep_find_examples
% %network_pair_info = load(fullfile(basedir,'helper_functions','network_pairs'));
% %network_pair_info = network_pair_info.network_pairs;
% %take this stuff directly from get_network_params()
% %look in get_network_params() for this NPjobs, 0 is the last one paired w/ 9..
% NPjobs = cat(1,[1:2:9],[2:2:9,0])'; %u-g-l-y
% num_net_types = size(NPjobs,1);
% network_pair_info = cell(num_net_types,1);
% for idx = 1:num_net_types
%     CP = num2cell(NPjobs(idx,:))';
%     CP = cellfun(@(x,y) {get_network_params(x,y)}, CP,repmat({struct()},2,1));
%     pair_table = cell2table(cell(2,numel(param_varnams)),'VariableNames',param_varnams);
%     pair_table.ItoE = cellfun(@(x) x.ItoE, CP,'UniformOutput',false);
%     pair_table.EtoI = cellfun(@(x) x.EtoI, CP,'UniformOutput',false);
%     pair_table.stim_A = cellfun(@(x) unique(x.trial_stimuli), CP,'UniformOutput',false);
%     pair_table.stim_B = cellfun(@(x) unique(x.trial_stimuli), CP,'UniformOutput',false);
%     pair_table.stim_targs = cellfun(@(x) x.stim_targs, CP,'UniformOutput',false);
%     pair_table.stim_targs = strrep(pair_table.stim_targs,'Estay','fast');
%     pair_table.stim_targs = strrep(pair_table.stim_targs,'Eswitch','slow');
%     network_pair_info{idx} = pair_table;
%     %     CP = cellfun(@(x) {x.ItoE,x.EtoI,unique(x.trial_stimuli),x.stim_targs},...
%     %         CP,'UniformOutput',false);
%     %     network_pair_info{idx} = cat(1,CP{:});
% end


fprintf('\n---loading simulation: %s\n',sim_name)

fprintf('\nHELLO need to check find_stay_durations() with verify argument')

%get results
fprintf('\nparfor disabled')
% num_workers = 24;
% c = parcluster('local');
% c.NumWorkers = num_workers;
% parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_stay_durations')})
% special_progress_tracker = fullfile(basedir,'SPT.txt');
% if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start

warning('\nfind_stay_durations() disabled, all non-undecided states')
file_data = cell(num_files,2);
%parfor idx = 1:num_files
for idx = 1:num_files
    %if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{1};
    
    %just get all of them, baseline test. Everything that's not undecided
    valid_states = ~strcmpi(state_durations(:,end),'undecided');
    %state.count recorded in second col 
    state_durations = state_durations(valid_states,2);
    state_durations = cat(1,state_durations{:});
    %convert to time
    state_durations = state_durations * timestep;
    %ratecheck estimates
    
%     state_durations = find_stay_durations(state_durations,curr_file.options);
%     switch pulse_stim
%         case 'yes' %just do this now while options is handy
%             state_durations = state_durations.samples;
%         case 'rem' %look at when IN the sample switch happened
%             state_durations = state_durations.decision_time;
%         case 'total_time'
%             state_durations = state_durations.duration;
%         case 'off'
%             state_durations = state_durations.duration;
%     end
    
    %store durations & parameters
    file_data(idx,:) = {state_durations,curr_file.options};
    %file_data{idx,1} = state_durations; %when this wasn't parfored
    %file_data{idx,2} = curr_file.options;
    

%     progress = worker_progress_tracker(special_progress_tracker);
%     if mod(progress,floor(num_files * .05)) == 0 %at half a percent
%         progress = (progress / num_files) * 100;
%         fprintf('----%.1f percent complete\n',progress);
%     end
end
%delete(gcp('nocreate'))
%delete(special_progress_tracker)

%search for jobs with identical parameters, collapse distributions
%get the randomized network parameters

switch opt.multiple_stimuli
    case 'yes'
        job_params = cellfun(@(x)...
            [x.ItoE, x.EtoI,x.trial_stimuli,find(strcmpi(x.stim_targs, stimtarg_vals))],...
            file_data(:,2),'UniformOutput',false); %matching "network_pair_info" format
    otherwise
        job_params = cellfun(@(x)...
            [x.ItoE, x.EtoI,unique(x.trial_stimuli),find(strcmpi(x.stim_targs, stimtarg_vals))],...
            file_data(:,2),'UniformOutput',false); %matching "network_pair_info" format
end

job_params = vertcat(job_params{:});
uniq_params = unique(job_params,'rows');
net_type = array2table(uniq_params,'VariableNames',param_varnams);
num_jobs = size(net_type,1);
fprintf('----------------------\n')
fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)

%collapse duplicate job parameters
result_data = cell(num_jobs,2);
Nruns = NaN(num_jobs,1); %record the number of successful jobs..
for idx = 1:num_jobs
    %find all matching
    curr_file = ismember(job_params,table2array(net_type(idx,:)),'rows');
    Nruns(idx) = sum(curr_file);
    explain_params = net_type(idx,:);
    explain_params.stim_targs = stimtarg_vals{explain_params.stim_targs};
    fprintf('\n---parameter set\n');disp(explain_params);fprintf('n files = %i\n',Nruns(idx))
    %collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    fprintf('------n states = %i\n',numel(result_data{idx,1}))
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
end

%this is stupid & obviously a hold-over from something I didn't implement well in the first place... 
net_type.stim_targs = cellfun(@(x) stimtarg_labels{x},...
    num2cell(net_type.stim_targs),'UniformOutput',false);

% %just grab some simple stats here
% mu_duration = cellfun(@(x) mean(x),result_data(:,1));
% med_duration = cellfun(@(x) median(x),result_data(:,1));
% logmu_dur = cellfun(@(x) mean(log(x)),result_data(:,1));
% %get the network parameters
% EtoE = cellfun(@(x)  x.EtoE,result_data(:,2));
% ItoE = cellfun(@(x)  x.ItoE,result_data(:,2));
% EtoI = cellfun(@(x)  x.EtoI,result_data(:,2));
% stim_value = cellfun(@(x)  unique(x.trial_stimuli),result_data(:,2));

% net_type = num2cell(num2cell(uniq_params),2);
% net_type = cellfun(@(x)  [x(1:3),stimtarg_vals{x{4}}],net_type,'UniformOutput',false);
% net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],net_type,'UniformOutput',false);
% net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],net_type,'UniformOutput',false);

keyboard


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
        Zlabel = sprintf('state duration (%s)',unit_measure);
    case 'med'
        Zlabel =  sprintf('median duration (%s)',unit_measure);
    case 'logmu'
        Zlabel = sprintf('log(%s) state duration',unit_measure);
        fig_fn = [fig_fn,'_log'];
end

if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;
alph = .5;

%for indexing the result paramters 
IDvars = param_varnams(~ismember(param_varnams,'stim_B')); %stim B not particular to network type

plt_idx = 0;
for idx = 1:num_net_types
    
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
        curr_data = table2cell(curr_net_info(j,IDvars));
        curr_data = cellfun(@(x) isequal(x,curr_data),... %ugly indexing & transform here..
            num2cell(table2cell(net_type(:,IDvars)),2));
        curr_data = result_data(curr_data,1);
        curr_data = cat(1,curr_data{:});
        
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = curr_data(curr_data ~= 0);%inf errors
                    curr_data = log(curr_data);
            end
            histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
            Ylab = 'freq';
        end
        %take control data
        BLinfo = curr_net_info(j,IDvars);
        BLinfo.stim_A = 0; BLinfo.stim_targs = 'baseline'; %0 stim & baseline targets 
        curr_data = cellfun(@(x) isequal(x,table2cell(BLinfo)),... 
            num2cell(table2cell(net_type(:,IDvars)),2));
        curr_data = result_data(curr_data,1);
        curr_data = cat(1,curr_data{:});
        
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
        switch pulse_stim
            case 'off'
                legend(sprintf('%.0fHz',curr_net_info.stim_A{j}),'location','best')
            otherwise
                legend(sprintf('%.0fHz %s',curr_net_info.stim_A{j},Pdur),'location','best')
        end
        
        if plt_idx == 9 || plt_idx == 10
            xlabel(Zlabel)
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info.stim_targs{j}),'Fontsize',14)
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
need_more = net_type(need_more,:);
current_count = num_states(num_states < 10000);
for idx = 1:num_types
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        curr_data = table2cell(curr_net_info(j,IDvars));
        curr_data = cellfun(@(x) isequal(x,curr_data),... %ugly indexing & transform here..
            num2cell(table2cell(need_more(:,IDvars)),2));
        if sum(curr_data) > 0
            curr_data = find(curr_data);
            fprintf('\nnetwork #%i %s has < 10k states',idx,curr_net_info.stim_targs{j})
            fprintf('\n---parameter sets:\n')
            for h = 1:numel(curr_data)
                disp(need_more(curr_data(h),:))
                fprintf('--- n = %i\n',current_count(curr_data(h)))
            end
            fprintf('\n------------------------\n')
        end
        
        BLinfo = curr_net_info(j,IDvars);
        BLinfo.stim_A = 0; BLinfo.stim_targs = 'baseline'; %0 stim & baseline targets
        curr_data = cellfun(@(x) isequal(x,table2cell(BLinfo)),... %ugly indexing & transform here..
            num2cell(table2cell(need_more(:,IDvars)),2));
        if sum(curr_data) > 0
            curr_data = find(curr_data);
            fprintf('\nnetwork #%i %s has < 10k states',idx,curr_net_info.stim_targs{j})
            fprintf('\n---parameter sets:\n')
            for h = 1:numel(curr_data)
                disp(need_more(curr_data(h),:))
                fprintf('--- n = %i\n',current_count(curr_data(h)))
            end
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


