clear
clc
format compact
hold off;close all


%opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
%opt.pulse_stim = 'rem'; %'yes' | 'total_time' | 'rem' whether to treat durations as samples (rem = time during sample)

Snames = {'dep_pref'};

figdir = 'figures_dep_pref'; %figdir = 'network_spiking_pulse_figures';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))


opt = struct();
opt.print_anything = 'yes'; %'yes' | 'no';
opt.load_or_calc = 'calc'; %load/calc (auto save with calc
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.rules = 'strict'; %strict|relaxed
opt.Tmax = NaN; %no transition limit
opt.pulse_stim = 'no'; % 'no' |'yes' | 'total_time' | 'rem' whether to treat durations as samples (rem = time during sample)
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'med';  %'mu' | 'med' | 'logmu'
make_my_figs(basedir,Snames{1},figdir,opt)


opt.load_or_calc = 'calc'; %load/calc (auto save with calc
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.Tmax = .3; %300 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'
opt.Tmax = .3; %300 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'med';  %'mu' | 'med' | 'logmu'
opt.Tmax = .3; %300 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

%broke here

opt.load_or_calc = 'calc'; %load/calc (auto save with calc
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
opt.Tmax = NaN; %no transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'med';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
make_my_figs(basedir,Snames{1},figdir,opt)


opt.load_or_calc = 'calc'; %load/calc (auto save with calc
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
opt.Tmax = .3; %250 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
opt.Tmax = .3; %250 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)

opt.load_or_calc = 'load'; %load/calc (auto save with calc
opt.outcome_stat = 'med';  %'mu' | 'med' | 'logmu'
opt.rules = 'relaxed'; %strict|relaxed
opt.Tmax = .3; %250 ms transition limit
%run preference analysis
make_my_figs(basedir,Snames{1},figdir,opt)





%
% %equate X axes for all figs of the same type
% %Snames = Snames(1:end-1);
% figdir =  fullfile(basedir,'Results',figdir,'preference');
% savedir = fullfile(figdir,'Xsynced');
% if ~isdir(savedir),mkdir(savedir);end
% ftypes = {'decision_timing_log','decision_timing','total_time','total_time_log','total_samples'};
%
% keyboard
% for idx = 1:numel(ftypes)
%
%     Xls = NaN(numel(Snames),2);
%     for fidx = 1:numel(Snames)
%         fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
%         h = openfig(fn);
%         Xls(fidx,:) = h.CurrentAxes.XLim;
%         close all
%     end
%     xlims = [min(Xls(:,1)),max(Xls(:,2))];
%     %xlims = [0,.8];
%     for fidx = 1:numel(Snames)
%         fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
%         h = openfig(fn);
%         ax = findobj(h,'Type','Axes');
%         set(ax,'XLim',xlims)
%
%         ch = allchild(ax);
%         ch = cat(1,ch{:});
%         set(ch,'Normalization','probability')
%
%         set(ch,'BinWidth',.01)
%         Fdir = fullfile(savedir,ftypes{idx});
%         if ~isdir(Fdir),mkdir(Fdir);end
%         print(fullfile(Fdir,[Snames{fidx} '_' ftypes{idx}]),'-djpeg')
%         close all
%     end
%
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
summary_stats = 'no'; %summary_stats = opt.summary_stats;

%handle directory & filename labels
figdir = fullfile(basedir,'Results',figdir,'preference');
svdir = fullfile(figdir,'data');
fig_fn = [sim_name '_%s_%s'];
Tfig_FN = 'timing_info_%s_%s';
svFN = sprintf('file_data_%s',opt.rules);
switch outcome_stat
    case 'mu'
        %do the main directory
        stat_label = 'mean';
    case 'med'
        stat_label = 'median';
    case 'logmu'
        stat_label = 'logmean';
end
figdir = fullfile(figdir,stat_label);
fig_fn = sprintf(fig_fn,opt.rules,stat_label);
Tfig_FN = sprintf(Tfig_FN,opt.rules,stat_label);
if ~isdir(figdir),mkdir(figdir);end
if ~isdir(svdir),mkdir(svdir);end
if ~isnan(opt.Tmax)
    label_Tmax = sprintf('_Tmax%i',round(opt.Tmax*1e3));
    fig_fn = [fig_fn, label_Tmax];
    Tfig_FN = [Tfig_FN, label_Tmax];
    svFN = [svFN, label_Tmax];
end
%get data filenames
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_BL'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);
num_files = numel(output_fns);
%labels for parameters and stimuli
param_varnams = {'ItoE','EtoI','stim_A','stim_B','stim_targs'};
stim_labels = {'A','B'};
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'};
timestep = .25e-3; %timestep is in set_options() now

switch pulse_stim
    case 'no'
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        Pdur = strsplit(sim_name,'_');
        Pdur = Pdur{3};
        Pdur = strrep(Pdur,'P','');
        Pdur = sprintf('%ss',Pdur);
end

%look in get_network_params() for this NPjobs, 0 is the last one paired w/ 9..
NPjobs = cat(1,[1:2:9],[2:2:9,0])'; %u-g-l-y
num_net_types = size(NPjobs,1);
network_pair_info = cell(num_net_types,1);
for idx = 1:num_net_types
    CP = num2cell(NPjobs(idx,:))';
    CP = cellfun(@(x,y) {get_network_params(x,y)}, CP,repmat({struct()},2,1));
    pair_table = cell2table(cell(2,numel(param_varnams)),'VariableNames',param_varnams);
    pair_table.ItoE = cellfun(@(x) x.ItoE, CP,'UniformOutput',false);
    pair_table.EtoI = cellfun(@(x) x.EtoI, CP,'UniformOutput',false);
    pair_table.stim_A = cellfun(@(x) unique(x.trial_stimuli), CP,'UniformOutput',false);
    pair_table.stim_B = cellfun(@(x) unique(x.trial_stimuli), CP,'UniformOutput',false);
    pair_table.stim_targs = cellfun(@(x) x.stim_targs, CP,'UniformOutput',false);
    pair_table.stim_targs = strrep(pair_table.stim_targs,'Estay','fast');
    pair_table.stim_targs = strrep(pair_table.stim_targs,'Eswitch','slow');
    network_pair_info{idx} = pair_table;
    %     CP = cellfun(@(x) {x.ItoE,x.EtoI,unique(x.trial_stimuli),x.stim_targs},...
    %         CP,'UniformOutput',false);
    %     network_pair_info{idx} = cat(1,CP{:});
end


fprintf('\n---loading simulation: %s\n',sim_name)

%get results
switch opt.load_or_calc
    case 'calc'
        num_workers = 24;
        num_workers = num_workers / 2;
        c = parcluster('local');
        c.NumWorkers = num_workers;
        parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_pref_durations')})
        special_progress_tracker = fullfile(basedir,'SPT.txt');
        if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start
        file_data = cell(num_files,3);
        parfor idx = 1:num_files
            %if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
            curr_file = load(output_fns{idx});
            %get state durations
            state_durations = curr_file.sim_results;
            state_durations = state_durations{1};
            
            [state_durations,Tinfo] = ...
                find_pref_durations(state_durations,curr_file.options,...
                'rules',opt.rules,'Tmax',opt.Tmax);
            
            if ~isempty(state_durations)
                Tstim = startsWith(state_durations.state,'stim');
            else
                Tstim = false;
            end
            state_durations = state_durations(Tstim,{'duration','state'});
            
            %I'm not handling pulse/vs continuous stimulus right now,
            %find_pref_durations() doesn't cover that either at the moment
            %     switch pulse_stim
            %         case 'yes' %just do this now while options is handy
            %             state_durations = state_durations.samples;
            %         case 'rem' %look at when IN the sample switch happened
            %             state_durations = state_durations.decision_time;
            %         case 'total_time'
            %             state_durations = state_durations.duration;
            %         case 'no'
            %             state_durations = state_durations.duration;
            %     end
            
            
            %store durations & parameters
            file_data(idx,:) = {state_durations,curr_file.options,Tinfo};
            
            progress = worker_progress_tracker(special_progress_tracker);
            if mod(progress,floor(num_files * .05)) == 0 %at half a percent
                progress = (progress / num_files) * 100;
                fprintf('----%.1f percent complete\n',progress);
            end
        end
        delete(gcp('nocreate'))
        delete(special_progress_tracker)
       
        %search for jobs with identical parameters, collapse distributions
        %get the randomized network parameters
        job_params = cellfun(@(x)...
            [x.ItoE, x.EtoI,x.trial_stimuli,find(strcmpi(x.stim_targs, stimtarg_vals))],...
            file_data(:,2),'UniformOutput',false); %matching "network_pair_info" format
        job_params = vertcat(job_params{:});
        uniq_params = unique(job_params,'rows');
        net_type = array2table(uniq_params,'VariableNames',param_varnams);
        num_jobs = size(net_type,1);
        fprintf('----------------------\n')
        fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)
        
        %collapse duplicate job parameters
        result_data = cell(num_jobs,2);
        timing_info = cell(num_jobs,1);
        Nruns = NaN(num_jobs,1); %record the number of successful jobs..
        for idx = 1:num_jobs
            %find all matching
            curr_file = ismember(job_params,table2array(net_type(idx,:)),'rows');
            Nruns(idx) = sum(curr_file);
            explain_params = net_type(idx,:);
            explain_params.stim_targs = stimtarg_vals{explain_params.stim_targs};
            fprintf('\n---parameter set\n');disp(explain_params);fprintf('n files = %i\n',Nruns(idx))
            %collapse & reallocate
            result_data{idx,1} = cat(1,file_data{curr_file,1});
            fprintf('------n states = %i\n',numel(result_data{idx,1}))
            %just grab the first options file... that shouldn't matter here
            result_data{idx,2} = file_data{find(curr_file,1),2};
            timing_info{idx} = cat(1,file_data{curr_file,3});
        end

        %save filedata
        save(fullfile(svdir,svFN),'result_data','timing_info','net_type','-v7.3')
        %clear out big clunky cell array 
        clear file_data 

    case 'load'
        %just load precomputed vars 
        result_data = load(fullfile(svdir,svFN));
        timing_info = result_data.timing_info;
        net_type = result_data.net_type;
        result_data = result_data.result_data;
end


%this is stupid & obviously a hold-over from something I didn't implement well in the first place...
net_type.stim_targs = cellfun(@(x) stimtarg_labels{x},...
    num2cell(net_type.stim_targs),'UniformOutput',false);

switch pulse_stim
    case 'yes'
        unit_measure = 'samples';
        %fig_fn = sprintf(fig_fn,'total_samples');
    case 'rem'
        unit_measure = 's - onset';
        %fig_fn = sprintf(fig_fn,'decision_timing');
    otherwise
        unit_measure = 's'; %like "seconds" not samples
        %fig_fn = sprintf(fig_fn,'total_time');
end

switch outcome_stat
    case 'mu'
        Ylabel = sprintf('state duration (%s)',unit_measure);
    case 'med'
        Ylabel =  sprintf('median duration (%s)',unit_measure);
    case 'logmu'
        Ylabel = sprintf('log(%s) state duration',unit_measure);
end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;
alph = .5;
fz = 10;

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
        alt_stim = net_type.stim_B(curr_data);
        curr_data = result_data(curr_data,1);
        
        if ~isempty(curr_data) %skip plot if no data...
            
            %make sure this is all sorted correctly
            [alt_stim,Bi] = sort(alt_stim);
            curr_data = curr_data(Bi);
            
            %split into cells for A & B (nested cellfuns are confusing but you can always verify...
            stim_data = cellfun(@(x) cellfun(@(y) strcmpi(sprintf('stim_%s',y),x.state),stim_labels,...
                'UniformOutput',false),curr_data,'UniformOutput',false);
            stim_data = cellfun(@(x,y)  cellfun(@(z) x(z,:),y,'UniformOutput',false),...
                curr_data,stim_data,'UniformOutput',false);
            stim_data = cat(1,stim_data{:});
            %cellfun(@(x) unique(x.state),stim_data) if you want to check again...
            switch outcome_stat
                case 'mu'
                    statfunc = @mean;
                case 'med'
                    statfunc = @median;
                case 'logmu'
                    %protect against inf errors too
                    statfunc = @(x) mean(log(x(x~=0)));
            end
            stim_stats = cellfun(@(x) varfun(statfunc,x(:,'duration'),'OutputFormat','uniform'),stim_data);
            stim_data = array2table([stim_stats,alt_stim],'VariableNames',[stim_labels,'Hz']);
            base_stim = curr_net_info.stim_A{j};
            
            Xvals = base_stim ./ stim_data.Hz;
            plot(Xvals,stim_data.A,'LineWidth',3)
            hold on
            plot(Xvals,stim_data.B,'LineWidth',3)
            %xlim([min(Xvals),max(Xvals)])
            
            legend_labs = {sprintf('%s: %.1fHz',stim_labels{1},base_stim),...
                sprintf('%s: varied',stim_labels{2})};
            legend(legend_labs,'Location','best','Orientation','horizontal','Box','off','Fontsize',fz)
            
        end
        
        hold off
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('A Hz / B Hz','Fontsize',fz,'FontWeight','bold')
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info.stim_targs{j}),'Fontsize',fz,'FontWeight','bold')
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('network #%i\n%s',idx,Ylabel),'Fontsize',fz,'FontWeight','bold')
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

figure()

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
        alt_stim = net_type.stim_B(curr_data);
        curr_data = timing_info(curr_data,1);
        
        if ~isempty(curr_data) %skip plot if no data...
            
            %make sure this is all sorted correctly
            [alt_stim,Bi] = sort(alt_stim);
            curr_data = curr_data(Bi);
            curr_data = cat(1,curr_data{:});
            
            %do log transform if logmu requested...
            switch outcome_stat
                case 'logmu'
                    %protect against inf errors too
                    statfunc = @(x) log(x(x~=0));
                    curr_data.u1 = statfunc(curr_data.u1);
                    curr_data.leave = statfunc(curr_data.leave);
                    curr_data.u2 = statfunc(curr_data.u2);
                    curr_data.total = statfunc(curr_data.total);
            end
            
            histogram(curr_data.total,'FaceColor',[.5,.5,.5],'EdgeColor',[.5,.5,.5],'FaceAlpha',alph);
            hold on
            histogram(curr_data.u1,'FaceAlpha',alph);
            histogram(curr_data.leave,'FaceAlpha',alph);
            histogram(curr_data.u2,'FaceAlpha',alph);
            
            
            if idx == 1 && j == 1
                axP = get(gca,'Position');
                lp = legend({'total','u1','leave','u2'},'FontSize',fz,...
                    'Location','northoutside','Box','off','Orientation','horizontal');
                set(gca, 'Position', axP)
            end
            lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];
        end
        
        hold off
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('A Hz / B Hz','Fontsize',fz,'FontWeight','bold')
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info.stim_targs{j}),'Fontsize',fz,'FontWeight','bold')
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('network #%i\n%s',idx,Ylabel),'Fontsize',fz,'FontWeight','bold')
        end
        
    end
end
orient tall
linkaxes(h,'x')
axis tight

switch print_anything
    case 'yes'
        print(fullfile(figdir,Tfig_FN),'-djpeg')
        savefig(fullfile(figdir,Tfig_FN))
end





%info about simulation
num_states = cellfun(@(x) numel(x{1}),num2cell(result_data,2));
need_more = num_states < 10000;
need_more = net_type(need_more,:);
current_count = num_states(num_states < 10000);
for idx = 1:num_net_types
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
        fprintf('%s:\n',Ylabel)
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


