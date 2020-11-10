clear
clc
format compact
hold off;close all
%note: if you ever wana make this nice... make a single function for
%creating a paramater ID cell array (or table) from an array of options structure.
%this function would be used for the template nets from get_network_params() and result files.
%Put this same function in network_spiking_results, etc. Also use cellfun(@(x) isequal(x,table2cell(curr_net_info(j,:))),Psets)

opt = struct();
opt.min_obs = 175; %min # of observations (states)
opt.print_anything = 'yes'; %'yes' | 'no';
opt.valid_states = 'stay'; %'stay' | 'all'; undecided is always invalid, 'all' gives stay & leave
opt.outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)
opt.parfor_load = 'on'; %on/off, must also (un)comment the actual for... line
opt.params2match = {'conn','stim'}; %!!!IMPORTANT!!! specify how results are matched to network types
%this can be at most {'conn','stim'}. That specifies matching on connection strengths, stimulus values


Snames = {'nets_D2t-slower_pref'};
figdir = cellfun(@(x) sprintf('figures_%s',x),Snames,'UniformOutput',false);


basedir = '~/Desktop/work/ACClab/rotation/project';
addpath(fullfile(basedir,'helper_functions'))


make_my_figs(basedir,Snames{1},figdir{1},opt);
return

for idx = 1:numel(Snames)
    opt.outcome_stat = 'mu';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt);
    opt.outcome_stat = 'logmu';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
    opt.outcome_stat = 'med';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
    opt.outcome_stat = 'logmed';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
end


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
opt.multiple_stimuli = 'yes'; %always yes...
outcome_stat = opt.outcome_stat;
pulse_stim = opt.pulse_stim;
print_anything = opt.print_anything; %'yes' | 'no';
summary_stats = 'no'; %summary_stats = opt.summary_stats;
params2match = opt.params2match;
figdir = fullfile(basedir,'Results',figdir,'durations');
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_baseline'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);
%checking for previously saved data
svdir = fullfile(figdir,'data');if ~isdir(svdir),mkdir(svdir);end
svFN = [sim_name '_%s.mat'];
switch pulse_stim
    case 'yes'
        svFN = sprintf(svFN,'total_samples');
    case 'rem'
        svFN = sprintf(svFN,'decision_timing');
    otherwise
        svFN = sprintf(svFN,'total_time');
end
if exist(fullfile(svdir,svFN)) > 0,load_summary = true;else,load_summary = false;end

switch opt.multiple_stimuli
    case 'yes'
        param_varnams = {'ItoE','EtoI','stim_A','stim_B','targ_cells'};
        %line below may be important
        %IDvars = param_varnams(~ismember(param_varnams,'stim_B')); %stim B not particular to network type
    case 'no'
        param_varnams = {'ItoE','EtoI','stim','targ_cells'};
end
%for indexing the result paramters
IDvars = [];
if sum(strcmp('conn',params2match)) > 0,IDvars = {'ItoE','EtoI'};end
if sum(strcmp('stim',params2match)) > 0,IDvars = [IDvars,param_varnams(startsWith(param_varnams,'stim'))];end

IDvars = IDvars(~ismember(IDvars,'stim_B')); %stim B not particular to network type

%get general options file from the first file
gen_options = load(output_fns{1});
gen_options = gen_options.options;
gen_options = rmfield(gen_options,{'stim_targs','trial_stimuli'});
timestep = gen_options.timestep;

switch pulse_stim
    case 'off'
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        error('get this from  options dude, was previously striped from sim name')
end

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'};

%info on the specific network parameters in this simulation
num_net_types = 10;
num_pairs = 5;
pair_inds = num2cell(reshape(1:num_net_types,[],num_pairs)); %just gives a cell array for pair indicies
network_pair_info = cell(num_pairs,1);
for idx = 1:num_pairs
    
    curr_params = cellfun(@(x) get_network_params(x,gen_options),pair_inds(:,idx),'UniformOutput',false);
    switch opt.multiple_stimuli
        case 'yes'
            curr_params = cellfun(@(x)...
                [x.ItoE, x.EtoI,num2cell(x.trial_stimuli{:}),x.stim_targs],...
                curr_params,'UniformOutput',false); %matching "network_pair_info" format
        otherwise
            curr_params = cellfun(@(x)...
                {x.ItoE, x.EtoI,unique(x.trial_stimuli{:}), x.stim_targs},...
                curr_params,'UniformOutput',false); %matching "network_pair_info" format
    end
    curr_params = cat(1,curr_params{:});
    T = cell2table(curr_params,'VariableNames',param_varnams);
    curr_types = T.targ_cells;
    curr_types = strrep(curr_types,'Estay','fast'); curr_types = strrep(curr_types,'Eswitch','slow');
    T.Properties.RowNames = curr_types;
    network_pair_info{idx} = T;
end

warning('find_pref_durations() will be needed for undecided states')


fprintf('\n---loading simulation: %s\n',sim_name)

if ~load_summary
    %get results
    switch opt.parfor_load
        case 'off'
            fprintf('\nparfor disabled\n')
        case 'on'
            num_workers = 24;
            c = parcluster('local');
            c.NumWorkers = num_workers;
            parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_stay_durations')})
            special_progress_tracker = fullfile(basedir,'SPT.txt');
            if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start
    end
    switch opt.valid_states %select states for analysis
        case 'all'
            warning('\nfind_stay_durations() disabled, all non-undecided states')
    end
    
    %warning('\ntrimming all data from events starting T < 20 s')
    
    file_data = cell(num_files,2);
    parfor idx = 1:num_files
        %for idx = 1:num_files
        switch opt.parfor_load
            case 'off'
                if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
        end
        curr_file = load(output_fns{idx});
        %get state durations
        state_durations = curr_file.sim_results;
        state_durations = state_durations{1};
        switch opt.valid_states %select states for analysis
            case 'all'
                %just get all of them, baseline test. Everything that's not undecided
                keep_states = ~strcmpi(state_durations(:,end),'undecided');
                %state.count recorded in second col
                state_durations = state_durations(keep_states,2);
                state_durations = cat(1,state_durations{:});
                %convert to time
                state_durations = state_durations * timestep;
            case 'stay'
                [state_durations,Sinfo] = find_stay_durations(state_durations,curr_file.options,'verify');
                
                %                 %limiting events to T > 20 s
                %                 too_early = cell2mat(Sinfo.event_time) - state_durations.duration; %start time
                %                 too_early = too_early < 20; %20 second limit
                %                 state_durations = state_durations(~too_early,:);
                %                 Sinfo = Sinfo(~too_early,:);
                
                switch pulse_stim
                    case 'yes' %just do this now while options is handy
                        state_durations = state_durations.samples;
                    case 'rem' %look at when IN the sample switch happened
                        state_durations = state_durations.decision_time;
                    case 'total_time'
                        state_durations = state_durations.duration;
                    case 'off'
                        state_durations = state_durations.duration;
                end
                state_durations = [array2table(state_durations,'VariableNames',{'data'}),Sinfo(:,'state')];
        end
        
        %store durations & parameters
        file_data(idx,:) = {state_durations,curr_file.options};
        
        switch opt.parfor_load
            case 'on'
                progress = worker_progress_tracker(special_progress_tracker);
                if mod(progress,floor(num_files * .05)) == 0 %at half a percent
                    progress = (progress / num_files) * 100;
                    fprintf('%s ---- %.1f percent complete\n',datestr(now,31),progress);
                end
        end
    end
    switch opt.parfor_load
        case 'on'
            delete(gcp('nocreate'))
            delete(special_progress_tracker)
    end
    
    fprintf('\nsaving data...\n')
    save(fullfile(svdir,svFN),'file_data')
elseif load_summary
    fprintf('\nloading saved summary data...\n')
    file_data = load(fullfile(svdir,svFN));
    file_data = file_data.file_data;
end

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
    explain_params.targ_cells = stimtarg_vals{explain_params.targ_cells};
    fprintf('\n---parameter set\n');disp(explain_params);fprintf('n files = %i\n',Nruns(idx))
    %collapse & reallocate
    this_data = file_data(curr_file,1); %so annoying...
    result_data{idx,1} = cat(1,this_data{:});
    fprintf('------n states = %i\n',size(result_data{idx,1},1))
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
end


Nobs = cellfun(@(x) numel(x(:,1)),result_data(:,1));
fprintf('\n\n:::: excluding sets wtih < %i observations\n',opt.min_obs)
min_obs = Nobs >= opt.min_obs;
fprintf('\n      %i sets excluded (out of %i total)\n\n',sum(~min_obs),numel(min_obs))
result_data = result_data(min_obs,:);
net_type = net_type(min_obs,:);

%this is stupid & obviously a hold-over from something I didn't implement well in the first place...
net_type.targ_cells = cellfun(@(x) stimtarg_labels{x},...
    num2cell(net_type.targ_cells),'UniformOutput',false);


wtf = result_data(:,2);
wtf = cellfun(@(x) abs(x.trial_stimuli(2) - 1092.9) < .05 && strcmp(x.stim_targs,'Estay'),wtf);
if sum(wtf) > 0,error('I forget why this is in here-- probably dataset specific');end
result_data = result_data(~wtf,:);
net_type = net_type(~wtf,:);

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
        Zlabel = sprintf('sampling (%s)',unit_measure);
    case 'med'
        Zlabel =  sprintf('median sampling (%s)',unit_measure);
        fig_fn = [fig_fn,'_med'];
    case 'logmu'
        Zlabel = 'seconds (log scale)';
        fig_fn = [fig_fn,'_log'];
    case 'logmed'
        Zlabel = sprintf('med. log_{10}(%s) sampling',unit_measure);
        fig_fn = [fig_fn,'_med_log'];
end

Ylab = 'p(x)';%Ylab = 'freq';

figdir = fullfile(figdir,sprintf('Nmin_%i',opt.min_obs)); %include the min observation cutoff
if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
fz = 20;
mk_sz = 300; %for adding netword symbs
mk_ln = 1;

%::::::::::::::::::::::Figure 6, supplementary:::::::::::::::::::::::::::::
%printing this figure first b/c it includes all data. This sets up the data
%for subsequent plots. 

net_symbs = {'o','square','^','diamond','v'}; %closed for fast nets, open for slow nets
%symbols to match parsweep_find_examples.m. Since get_network_params() creates
%network_pair_info, the ordering below  will match the network_pair_info
%ordering in parsweep_find_examples.m.
SPdata = []; %for scatter plot data
h = [];
plt_idx = 0;
figure;orient tall
for idx = 1:num_pairs
    
    curr_net_info = network_pair_info{idx};
    for j = 2:-1:1 %match the fast, slow ordering in other figures...
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(ceil(num_net_types/2),2,plt_idx);
        hold on
        
        %find the right results for network set-up
        net_ind = curr_net_info{j,IDvars};
        net_ind = ismember(net_type{:,IDvars},net_ind,'rows');
        curr_data = result_data(net_ind,1);
        
        curr_type = curr_net_info.Row{j};
        switch outcome_stat
            case 'mu'
                statfunc = @mean;
            case 'med'
                statfunc = @median;
            case 'logmu'
                %protect against inf errors too
                statfunc = @(x) mean(log10(x(x~=0)));
            case 'logmed'
                statfunc = @(x) median(log10(x(x~=0)));
        end
        
        curr_data = cellfun(@(x) varfun(statfunc,x,'InputVariables','data',...
            'GroupingVariables','state') ,curr_data,'UniformOutput',false);
        curr_data = cellfun(@(x) x(:,[1,size(x,2)]),curr_data,'UniformOutput',false);
        curr_data = cellfun(@(x) array2table(x{:,size(x,2)}','Variablenames',strrep(x.state,'stim','data')),...
            curr_data,'UniformOutput',false); %turn into 1 x 2 table w/ stim A/B as varnames
        curr_data = cat(1,curr_data{:});
        %now take the net info as well, so it's easy
        curr_data = [net_type(net_ind,:),curr_data];
        curr_data = sortrows(curr_data,'stim_B'); %sort by alternate stim strength
        
        base_stim = curr_net_info.stim_A(j);
        Xvals = curr_data.stim_B ./ base_stim;
        
        switch curr_type
            case 'slow'
                col = matorange;
                nm = scatter(NaN,NaN,mk_sz,'black',net_symbs{idx},...
                    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'LineWidth',mk_ln);
                legend(nm,'network','Location','northeast','AutoUpdate','off')
                tit = 'slow (repel) networks';
            case 'fast'
                col = matblue;
                nm = scatter(NaN,NaN,mk_sz,'black',net_symbs{idx},'filled',...
                    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
                legend(nm,'network','Location','northwest','AutoUpdate','off')
                tit = 'fast (entice) networks';
        end
        
        
        plot(Xvals,curr_data.data_A,'-o','LineWidth',2,'Color',col)
        hold on
        plot(Xvals,curr_data.data_B,':o','LineWidth',2,'Color',col)
        
        %xlim([min(Xvals),max(Xvals)])
       
        hold off
        axis tight
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('B / Hz','FontWeight','bold')
        end
        if plt_idx == 1 || plt_idx == 2
            title(tit,'FontWeight','bold','Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(Zlabel,'FontWeight','bold')
        end
        
        %organize data for scatter plo
        curr_data.net_index = repmat(plt_idx,size(curr_data,1),1); %index ID for scatter plot
        SPdata = [SPdata;curr_data];
    end
end


switch outcome_stat
    case {'logmu','logmed'}
        linkaxes(h,'y')
        linkaxes(h(1:2:end),'x');linkaxes(h(2:2:end),'x')
        
        for idx = 1:numel(h)
            axes(h(idx))
            Ytick = get(gca,'YTick');
            %Ytick = linspace(Ytick(1),Ytick(end),5);
            %set(gca,'YTick',Ytick)
            Ytick = 10.^Ytick;
            Ytick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Ytick),'UniformOutput',false);
            Ytick = strrep(Ytick,'.0s','s');
            Ytick = strrep(Ytick,'0.','.');
            Ytick = strrep(Ytick,'s',''); %went from "sampling" label to just seconds...
            set(gca,'YTickLabel',Ytick);
        end
        
        
    case {'mu','med'}
        linkaxes(h(1:2:end),'x');linkaxes(h(2:2:end),'x')
        %axis tight
end



switch print_anything
    case 'yes'
        set(gcf,'Renderer','painters')
        %my code typically saves figures in specific results directories, w/ particular filenames
        print(fullfile(figdir,fig_fn),'-djpeg','-r600')
        savefig(fullfile(figdir,fig_fn))
        %also save one here since it's a figure for the paper
        print('fig6-full_data','-djpeg','-r600')
end


%::::::::::::::::::::::Figure 5::::::::::::::::::::::::::::::::::::::::::::

close all;figure;orient portrait
%scatter plot
alph = .75;Msz = 75;
SPdata.X = SPdata.stim_B ./ SPdata.stim_A;
SPcells = {'fast','slow'};
leg_labs = {'fast (entice) network','slow (repel) network'};

for idx = 1:numel(SPcells)
    hold on
    curr_data = strcmp(SPdata.targ_cells,SPcells{idx});
    curr_data = SPdata(curr_data,:);
    curr_nets = unique(curr_data.net_index);
    %for plt_idx = 1:numel(curr_nets)
    for plt_idx = 2
        this_net = curr_data.net_index == curr_nets(plt_idx);
        this_net = curr_data(this_net,:);
        switch SPcells{idx}
            case 'slow'
                col = matorange;
            case 'fast'
                col = matblue;
        end
        
        %scatter(this_net.data_B,this_net.data_A,Msz,col,'Marker',markers{plt_idx},'Linewidth',1.25)
        scatter(this_net.data_B,this_net.data_A,Msz,col,'filled','MarkerFaceAlpha',alph) %slightly different coloring/markers here
        
        axis tight
        
        switch outcome_stat
            case 'logmu'
                xlabel({'B - varied stimuli','seconds (log scale)'},'FontWeight','bold')
                ylabel({'A - constant stimuli','seconds (log scale)'},'FontWeight','bold')
            otherwise
                xlabel(sprintf(['B - varied [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
                ylabel(sprintf(['A - constant [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
        end
    end
end

%title(sprintf('Implicit Competition'),'FontWeight','bold','Fontsize',14)

set(gca,'FontSize',fz)
Xtick = get(gca,'XTick');
Xtick = linspace(Xtick(1),Xtick(end),5);
set(gca,'XTick',Xtick)
Xtick = 10.^Xtick;
Xtick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Xtick),'UniformOutput',false);
Xtick = strrep(Xtick,'.0s','s');
Xtick = strrep(Xtick,'0.','.');
Xtick = strrep(Xtick,'s',''); %went from "sampling" label to just seconds...
set(gca,'XTickLabel',Xtick);

Ytick = get(gca,'YTick');
Ytick = linspace(Ytick(1),Ytick(end),5);
set(gca,'YTick',Ytick)
Ytick = 10.^Ytick;
Ytick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Ytick),'UniformOutput',false);
Ytick = strrep(Ytick,'.0s','s');
Ytick = strrep(Ytick,'0.','.');
Ytick = strrep(Ytick,'s',''); %went from "sampling" label to just seconds...
set(gca,'YTickLabel',Ytick);

[~,hobj] = legend(leg_labs,'Box','off','Location','southwest','FontSize',fz);
ll = findobj(hobj,'type','patch');
set(ll,'MarkerSize',sqrt(Msz),'FaceAlpha',alph);

switch print_anything
    case 'yes'
        set(gcf,'Renderer','painters')
        print('fig5-preference_data','-djpeg','-r600') %schwartzupdate_fig
end

%::::::::::::::::::::::Figure 6::::::::::::::::::::::::::::::::::::::::::::

%This is for the fig 6 lineplots 
h = [];
plt_idx = 0;
figure;set(gcf,'Renderer','painters')
scr = get( groot,'Screensize');
fwid = 1000;
pos = [1,scr(3),fwid,fwid*.4];
set(gcf,'Position',pos);
for idx = 2
    
    curr_net_info = network_pair_info{idx};
    for j = 2:-1:1 %match the fast, slow ordering in other figures...
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(1,2,plt_idx);
        hold on
        
        %find the right results for network set-up
        net_ind = curr_net_info{j,IDvars};
        net_ind = ismember(net_type{:,IDvars},net_ind,'rows');
        curr_data = result_data(net_ind,1);
        
        curr_type = curr_net_info.Row{j};
        switch outcome_stat
            case 'mu'
                statfunc = @mean;
            case 'med'
                statfunc = @median;
            case 'logmu'
                %protect against inf errors too
                statfunc = @(x) mean(log10(x(x~=0)));
            case 'logmed'
                statfunc = @(x) median(log10(x(x~=0)));
        end
        
        curr_data = cellfun(@(x) varfun(statfunc,x,'InputVariables','data',...
            'GroupingVariables','state') ,curr_data,'UniformOutput',false);
        curr_data = cellfun(@(x) x(:,[1,size(x,2)]),curr_data,'UniformOutput',false);
        curr_data = cellfun(@(x) array2table(x{:,size(x,2)}','Variablenames',strrep(x.state,'stim','data')),...
            curr_data,'UniformOutput',false); %turn into 1 x 2 table w/ stim A/B as varnames
        curr_data = cat(1,curr_data{:});
        %now take the net info as well, so it's easy
        curr_data = [net_type(net_ind,:),curr_data];
        curr_data = sortrows(curr_data,'stim_B'); %sort by alternate stim strength
        
        base_stim = curr_net_info.stim_A(j);
        Xvals = curr_data.stim_B ./ base_stim;
        
        switch curr_type
            case 'slow'
                col = matorange;
                tit = 'slow (repel) network';
                
            case 'fast'
                col = matblue;
                tit = 'fast (entice) network';
                
        end
        
        plot(Xvals,curr_data.data_A,'LineWidth',3,'Color',col)
        hold on
        plot(Xvals,curr_data.data_B,'LineWidth',3,'Color',col,'LineStyle','--')
        %xlim([min(Xvals),max(Xvals)])
        %legend_labs = {sprintf('%s: %.0f Hz','A',base_stim),...
        %    sprintf('%s: varied','B')};
        %legend_labs = {'A - constant stimuli','B - varied stimuli'};
        
        %         sampling_change = curr_data{[1,size(curr_data,1)],{'data_A','data_B'}}; %beginning & end
        %         sampling_change = diff(sampling_change);
        %         %legend_labs = cellfun(@(x,y) [x '\newline\Deltay = ' sprintf('%.2f',y)],...
        %         %    legend_labs,num2cell(sampling_change),'UniformOutput',false);
        %         legend_labs = cellfun(@(x,y) [x ' (\Deltay = ' sprintf('%.2f)',y)],...
        %             legend_labs,num2cell(sampling_change),'UniformOutput',false);
        
        
        legend_labs = {'A - constant stimuli','B - varied stimuli'};
        legend(legend_labs,'Location','northwest','Box','off')
        
        
        axis tight
        
        switch curr_type
            case 'slow'
                xlim([0,2])
                legend(legend_labs,'Location','northeast','Box','off')
        end
        
        
        set(gca,'FontSize',fz)
        xlabel('B / A ','FontWeight','bold')
        switch outcome_stat
            case 'logmu'
                ylabel('seconds (log scale)','FontWeight','bold')
            otherwise
                ylabel(Zlabel,'FontWeight','bold')
        end
        title(tit,'FontWeight','bold')
        
    end
end

linkaxes(h,'y')
switch outcome_stat
    case 'logmu'
        for idx = 1:numel(h)
            axes(h(idx))
            Ytick = get(gca,'YTick');
            Ytick = linspace(Ytick(1),Ytick(end),5);
            set(gca,'YTick',Ytick)
            Ytick = 10.^Ytick;
            Ytick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Ytick),'UniformOutput',false);
            Ytick = strrep(Ytick,'.0s','s');
            Ytick = strrep(Ytick,'0.','.');
            Ytick = strrep(Ytick,'s',''); %went from "sampling" label to just seconds...
            set(gca,'YTickLabel',Ytick);
        end
end

switch print_anything
    case 'yes'
        set(gcf,'Renderer','painters')
        print('fig6-line_data','-djpeg','-r600')
end


%::::::::::::::::::::::Figure 5, supplementary:::::::::::::::::::::::::::::


close all;figure;orient portrait
mk_ln = 1.5; %mk_sz = 75
fz = 20;
for idx = 1:numel(SPcells)
    hold on
    curr_data = strcmp(SPdata.targ_cells,SPcells{idx});
    curr_data = SPdata(curr_data,:);
    curr_nets = unique(curr_data.net_index);
    for plt_idx = 1:numel(curr_nets)
        this_net = curr_data.net_index == curr_nets(plt_idx);
        this_net = curr_data(this_net,:);
        switch SPcells{idx}
            case 'slow'
                col = matorange;
                scatter(this_net.data_B,this_net.data_A,Msz,col,...
                    net_symbs{plt_idx},'MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph,'LineWidth',mk_ln);
                
            case 'fast'
                col = matblue;
                scatter(this_net.data_B,this_net.data_A,Msz,col,'filled',...
                    net_symbs{plt_idx},'MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph);

        end

        axis tight
        
        switch outcome_stat
            case 'logmu'
                xlabel({'B - varied stimuli','seconds (log scale)'},'FontWeight','bold')
                ylabel({'A - constant stimuli','seconds (log scale)'},'FontWeight','bold')
            otherwise
                xlabel(sprintf(['B - varied [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
                ylabel(sprintf(['A - constant [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
        end
    end
end

set(gca,'FontSize',fz)
Xtick = get(gca,'XTick');
Xtick = linspace(Xtick(1),Xtick(end),5);
set(gca,'XTick',Xtick)
Xtick = 10.^Xtick;
Xtick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Xtick),'UniformOutput',false);
Xtick = strrep(Xtick,'.0s','s');
Xtick = strrep(Xtick,'0.','.');
Xtick = strrep(Xtick,'s',''); %went from "sampling" label to just seconds...
set(gca,'XTickLabel',Xtick);

Ytick = get(gca,'YTick');
Ytick = linspace(Ytick(1),Ytick(end),5);
set(gca,'YTick',Ytick)
Ytick = 10.^Ytick;
Ytick = cellfun(@(x) sprintf('%.1fs',x),num2cell(Ytick),'UniformOutput',false);
Ytick = strrep(Ytick,'.0s','s');
Ytick = strrep(Ytick,'0.','.');
Ytick = strrep(Ytick,'s',''); %went from "sampling" label to just seconds...
set(gca,'YTickLabel',Ytick);

% %need something else here!
% [~,hobj] = legend(strcat(SPcells,' network'),'Box','off','Location','southwest','FontSize',20);
% ll = findobj(hobj,'type','patch');
% set(ll,'MarkerSize',sqrt(Msz),'FaceAlpha',alph);

% switch print_anything
%     case 'yes'
%         set(gcf,'Renderer','painters')
%         print('fig5-preference_data','-djpeg','-r600') %schwartzupdate_fig
% end

%this is for the nice legend with all the different colors & symbols------
lg_pos = legend(' ','Location','best'); 
if contains(lg_pos.Location,'east')
    make_lg = 'right';
else
    make_lg = 'left';
end
if contains(lg_pos.Location,'north')
    Y0 = .95;
else
    Y0 = .2;
end
delete(lg_pos)
switch make_lg
    case 'left'
        X0 = .025; move_pt = .02;
    case 'right'
        X0 = .975; move_pt = -.02;%was   Y0 = .5
end
lg_fz = fz - 2;
%leg_labs = strcat(SPcells,' networks');
leg_labs = strcat(leg_labs,'s');
get_ax_val = @(p,x) p*range(x)+min(x);
xlim(xlim + [-(range(xlim)*.015),(range(xlim)*.015)]); %extra room
ylim(ylim + [-(range(ylim)*.015),(range(ylim)*.015)]);
X = get_ax_val(X0,xlim);
text(X,get_ax_val(Y0,ylim),leg_labs{1},'Fontsize',lg_fz,'HorizontalAlignment',make_lg)
text(X,get_ax_val(Y0-.1,ylim),leg_labs{2},'Fontsize',lg_fz,'HorizontalAlignment',make_lg)
%X = get_ax_val(.975,xlim);
%text(X,get_ax_val(.5,ylim),SPcells{1},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment','right')
%text(X,get_ax_val(.4,ylim),SPcells{2},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment','right')
for plt_idx = 1:numel(curr_nets)
    
    %X = get_ax_val(.975-(plt_idx*.02),xlim);
    %scatter(X,get_ax_val(.45,ylim),Msz,Bcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    %scatter(X,get_ax_val(.35,ylim),Msz,Rcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    X = get_ax_val(X0+(plt_idx*move_pt),xlim);
    
    % %regular w/ colors
    % scatter(X,get_ax_val(Y0-.05,ylim),Msz,Bcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    % scatter(X,get_ax_val(Y0-.15,ylim),Msz,Rcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    
    %for symbols
    scatter(X,get_ax_val(Y0-.05,ylim),Msz,matblue,'filled',net_symbs{plt_idx},...
        'MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph);
    scatter(X,get_ax_val(Y0-.15,ylim),Msz,matorange,net_symbs{plt_idx},...
        'MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph,'LineWidth',mk_ln);
end

%above for nice legend------
switch print_anything
    case 'yes'
        set(gcf,'Renderer','painters')
        %my code typically saves figures in specific results directories, w/ particular filenames
        print(fullfile(figdir,[fig_fn '-scatter']),'-djpeg','-r600')
        savefig(fullfile(figdir,[fig_fn '-scatter']))
        %also save one here since it's a figure for the paper
        print('fig5-full_data','-djpeg','-r600')
end


%info about simulation
num_states = cellfun(@(x) size(x{1},1),num2cell(result_data,2));
need_more = num_states < 10000;
need_more = net_type(need_more,:);
current_count = num_states(num_states < 10000);
for idx = 1:num_pairs
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        curr_data = table2cell(curr_net_info(j,IDvars));
        curr_data = cellfun(@(x) isequal(x,curr_data),... %ugly indexing & transform here..
            num2cell(table2cell(need_more(:,IDvars)),2));
        if sum(curr_data) > 0
            curr_data = find(curr_data);
            fprintf('\nnetwork #%i %s has < 10k states',idx,curr_net_info.targ_cells{j})
            fprintf('\n---parameter sets:\n')
            for h = 1:numel(curr_data)
                disp(need_more(curr_data(h),:))
                fprintf('--- n = %i\n',current_count(curr_data(h)))
            end
            fprintf('\n------------------------\n')
        end
        
        BLinfo = curr_net_info(j,IDvars);
        BLinfo{:,startsWith(param_varnams,'stim')} = 0; BLinfo.targ_cells = 'baseline';
        curr_data = cellfun(@(x) isequal(x,table2cell(BLinfo)),... %ugly indexing & transform here..
            num2cell(table2cell(need_more),2));
        if sum(curr_data) > 0
            curr_data = find(curr_data);
            fprintf('\nnetwork #%i %s has < 10k states',idx,curr_net_info.targ_cells{j})
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
                        curr_data = log10(curr_data(curr_data~=0));
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
                    Xstim = cellfun(@(x) log10(x(x~=0)),Xstim,'UniformOutput',false);
            end
            [~,pval] = ttest2(Xstim{1},Xstim{2}); %regular old t-test
            fprintf('t-test p = %.3f\n',pval)
            [CI,H] = boot_mean_diffs(Xstim{1},Xstim{2},10000);
            fprintf('bootstrap test: %s\n',H)
            fprintf('bootstrap CI: %.2f, %.2f\n',CI)
        end
end
end


