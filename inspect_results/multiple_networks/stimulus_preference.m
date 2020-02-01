clear
clc
format compact
hold off;close all
%note: if you ever wana make this nice... make a single function for
%creating a paramater ID cell array (or table) from an array of options structure.
%this function would be used for the template nets from get_network_params() and result files.
%Put this same function in network_spiking_results, etc. Also use cellfun(@(x) isequal(x,table2cell(curr_net_info(j,:))),Psets)

opt = struct();
opt.print_anything = 'yes'; %'yes' | 'no';
opt.valid_states = 'stay'; %'stay' | 'all'; undecided is always invalid, 'all' gives stay & leave
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)
opt.parfor_load = 'on'; %on/off, must also (un)comment the actual for... line
opt.params2match = {'conn','stim'}; %!!!IMPORTANT!!! specify how results are matched to network types
%this can be at most {'conn','stim'}. That specifies matching on connection strengths, stimulus values


Snames = {'nets_D2t-slower_pref'};
figdir = cellfun(@(x) sprintf('figures_%s',x),Snames,'UniformOutput',false);


%basedir = '~/Desktop/ksander/rotation/project';
basedir = '~/Desktop/work/ACClab/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

for idx = 1:numel(Snames)
    opt.outcome_stat = 'mu';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
    opt.outcome_stat = 'logmu';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
    opt.outcome_stat = 'med';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
    opt.outcome_stat = 'logmed';
    make_my_figs(basedir,Snames{idx},figdir{idx},opt)
end
return

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

% %control, log total time @ 150 pulse
% opt.outcome_stat = 'mu'; opt.pulse_stim = 'off';
% make_my_figs(basedir,Snames{1},figdir,opt)
%
% %control, log total time @ 10 pulse
% opt.outcome_stat = 'logmu'; opt.pulse_stim = 'total_time';
% make_my_figs(basedir,Snames{1},figdir,opt)


%equate X axes for all figs of the same type
%Snames = Snames(1:end-1);
figdir = fullfile(basedir,'Results',figdir,'durations');
savedir = fullfile(figdir,'Xsynced');
if ~isdir(savedir),mkdir(savedir);end
ftypes = {'total_time','total_time_log'};
%ftypes = {'decision_timing_log','decision_timing','total_time','total_time_log','total_samples'};

invid_cols = {'yes','no'}; %equate x axes within column (e.g. slow vs fast)

for idx = 1:numel(ftypes)
    
    Xls = NaN(numel(Snames),2);
    switch invid_cols{idx}
        case 'yes'
            Xls = cat(3,Xls,Xls);
    end
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        switch invid_cols{idx}
            case 'yes'
                ax = num2cell(findall(h, 'type', 'axes'));
                %left column even, right column odd
                ax = cellfun(@(x) x.XLim,ax,'UniformOutput',false);
                ax = [ax(2:2:end),ax(1:2:end)]; %sort to columns
                ax = num2cell(ax,1);
                ax = cellfun(@cell2mat,ax,'UniformOutput',false);
                ax = cellfun(@(x) [min(x(:,1)),max(x(:,2))],ax,'UniformOutput',false);
                Xls(fidx,:,:) = cat(3,ax{:}); %put in 3rd D
            otherwise
                Xls(fidx,:) = h.CurrentAxes.XLim;
        end
        close all
    end
    switch invid_cols{idx}
        case 'yes'
            xlims = [min(Xls(:,1,:)),max(Xls(:,2,:))];
            xlims = squeeze(xlims); %now Xlims column is Xlims for a column
        otherwise
            xlims = [min(Xls(:,1)),max(Xls(:,2))];
    end
    %xlims = [0,.8];
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        ax = findobj(h,'Type','Axes');
        switch invid_cols{idx}
            case 'yes' %left column even, right column odd
                set(ax(2:2:end),'XLim',xlims(:,1))
                set(ax(1:2:end),'XLim',xlims(:,2))
            otherwise
                set(ax,'XLim',xlims)
        end
        ch = allchild(ax);
        ch = cat(1,ch{:});
        set(ch,'Normalization','probability')
        binwit = cellfun(@(x) x.XLim,num2cell(ax),'UniformOutput',false);
        binwit = cellfun(@range,binwit) ./ 50;
        binwit = num2cell(binwit);
        nodata = cellfun(@(x) isempty(x.Children),num2cell(ax));
        cellfun(@(x,y) set(x,'BinWidth',y),num2cell(ch),binwit(~nodata));   %set(ch,'BinWidth',.1)
        Fdir = fullfile(savedir,ftypes{idx});
        if ~isdir(Fdir),mkdir(Fdir);end
        print(fullfile(Fdir,[Snames{fidx} '_' ftypes{idx}]),'-djpeg')
        close all
    end
    
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
                [x.ItoE, x.EtoI,num2cell(x.trial_stimuli),x.stim_targs],...
                curr_params,'UniformOutput',false); %matching "network_pair_info" format
        otherwise
            curr_params = cellfun(@(x)...
                {x.ItoE, x.EtoI,unique(x.trial_stimuli), x.stim_targs},...
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
            warning('parfor set for bender: # workers = 4') 
            num_workers = 4; %24;
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

%this is stupid & obviously a hold-over from something I didn't implement well in the first place...
net_type.targ_cells = cellfun(@(x) stimtarg_labels{x},...
    num2cell(net_type.targ_cells),'UniformOutput',false);


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
        Zlabel = sprintf('log_{10}(%s) sampling',unit_measure);
        fig_fn = [fig_fn,'_log'];
    case 'logmed'
        Zlabel = sprintf('med. log_{10}(%s) sampling',unit_measure);
        fig_fn = [fig_fn,'_med_log'];
end

Ylab = 'p(x)';%Ylab = 'freq';

if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;

SPdata = []; %for scatter plot data
h = [];
plt_idx = 0;
for idx = 1:num_pairs
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(ceil(num_net_types/2),2,plt_idx);
        hold on
                
        %find the right results for network set-up
        net_ind = curr_net_info{j,IDvars};
        net_ind = ismember(net_type{:,IDvars},net_ind,'rows');
        curr_data = result_data(net_ind,1);
        
        %         %this is just to save data for inspection...
        %         data_table = cell(size(curr_data));
        %         for CC = 1:numel(curr_data)
        %             CCC = curr_data{CC};
        %             CCC = cell2table(cellfun(@(x) CCC.data(ismember(CCC.state,x)) ,unique(CCC.state),'UniformOutput',false)',...
        %                 'VariableNames',strrep(unique(CCC.state)','stim','data'));
        %             data_table{CC} = CCC;
        %         end
        %         data_table = [net_type(net_ind,:),cat(1,data_table{:})];
        %         save(fullfile(svdir,sprintf('%s_net%i',curr_net_info.Row{j},idx)),'data_table')
        %         %end inspection saving code...
        
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
        
        plot(Xvals,curr_data.data_A,'LineWidth',3)
        hold on
        plot(Xvals,curr_data.data_B,'LineWidth',3)
        %xlim([min(Xvals),max(Xvals)])
        legend_labs = {sprintf('%s: %.0f Hz','A',base_stim),...
            sprintf('%s: varied','B')};
                
        sampling_change = curr_data{[1,size(curr_data,1)],{'data_A','data_B'}}; %beginning & end
        sampling_change = diff(sampling_change);
        %legend_labs = cellfun(@(x,y) [x '\newline\Deltay = ' sprintf('%.2f',y)],...
        %    legend_labs,num2cell(sampling_change),'UniformOutput',false);
        legend_labs = cellfun(@(x,y) [x ' (\Deltay = ' sprintf('%.2f)',y)],...
            legend_labs,num2cell(sampling_change),'UniformOutput',false);
        legend(legend_labs,'Location','best','Box','off')
        
        
        hold off
        axis tight
        %legend(sprintf('\\mu = %.1f',mean(curr_data)),'location','best')
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('B Hz / A Hz','FontWeight','bold')
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info.Row{j}),'FontWeight','bold','Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('network #%i\n%s',idx,Zlabel),'FontWeight','bold')
        end
        
        %organize data for scatter plot
        curr_data.net_index = repmat(plt_idx,size(curr_data,1),1); %index ID for scatter plot
        SPdata = [SPdata;curr_data];
    end
end

orient tall
switch outcome_stat
    case {'logmu','logmed'}
        linkaxes(h,'y')
    case {'mu','med'}
        linkaxes(h(1:2:end),'x');linkaxes(h(2:2:end),'x')
        axis tight
end

switch print_anything
    case 'yes'
        print(fullfile(figdir,fig_fn),'-djpeg')
        savefig(fullfile(figdir,fig_fn))
end
close all;figure;orient portrait
%scatter plot 
alph = .75;Msz = 75;
Bcol = colormap('winter');Rcol = colormap('autumn');
SPdata.X = SPdata.stim_B ./ SPdata.stim_A;
SPcells = {'slow','fast'};
%markers = {'o','square','diamond','pentagram','hexagram'};
markers = {'x','+','^','v','d'};
for idx = 1:numel(SPcells)
    %subplot(numel(SPcells),1,idx);hold on
    hold on
    curr_data = strcmp(SPdata.targ_cells,SPcells{idx});
    curr_data = SPdata(curr_data,:);
    curr_nets = unique(curr_data.net_index);
    for plt_idx = 1:numel(curr_nets)
        this_net = curr_data.net_index == curr_nets(plt_idx);
        this_net = curr_data(this_net,:);
        col_idx = floor(size(Bcol,1)./numel(curr_nets)).*(plt_idx-1) + 1; %color index
        switch SPcells{idx}
            case 'slow'
                col = Bcol(col_idx,:);
            case 'fast'
                col = Rcol(col_idx,:);
        end
        
        %scatter(this_net.X,this_net.data_A,Msz,col,'filled','MarkerFaceAlpha',alph)
        %scatter(this_net.data_B,this_net.data_A,Msz,col,'filled','MarkerFaceAlpha',alph)
        %scatter(this_net.data_B,this_net.data_A,Msz,col,'filled','MarkerFaceAlpha',alph,'Marker',markers{plt_idx})
        scatter(this_net.data_B,this_net.data_A,Msz,col,'Marker',markers{plt_idx})
        axis tight
        title(sprintf('Implicit Competition'),...
            'FontWeight','bold','Fontsize',14)
        %ylabel(Zlabel,'FontWeight','bold')
        
        %xlabel('proportion of alternative stimulus','FontWeight','bold')
        
        xlabel(sprintf(['B - varied [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
        ylabel(sprintf(['A - constant [' strrep(Zlabel,' sampling','') ']']),'FontWeight','bold')
    end
end
make_lg = 'left';
switch make_lg
    case 'left'
        Y0 = .2; X0 = .025; move_pt = .02;
    case 'right'
        Y0 = .2; X0 = .975; move_pt = -.02;%was   Y0 = .5
end
lg_fz = 16; get_ax_val = @(p,x) p*range(x)+min(x);
xlim(xlim + [-(range(xlim)*.015),(range(xlim)*.015)]); %extra room
ylim(ylim + [-(range(ylim)*.015),(range(ylim)*.015)]);
X = get_ax_val(X0,xlim);
text(X,get_ax_val(Y0,ylim),SPcells{1},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment',make_lg) 
text(X,get_ax_val(Y0-.1,ylim),SPcells{2},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment',make_lg) 
%X = get_ax_val(.975,xlim);
%text(X,get_ax_val(.5,ylim),SPcells{1},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment','right') 
%text(X,get_ax_val(.4,ylim),SPcells{2},'Fontsize',lg_fz,'FontWeight','bold','HorizontalAlignment','right') 
for plt_idx = 1:numel(curr_nets)
    col_idx = floor(size(Bcol,1)./numel(curr_nets)).*(plt_idx-1) + 1; %color index
    %X = get_ax_val(.975-(plt_idx*.02),xlim);
    %scatter(X,get_ax_val(.45,ylim),Msz,Bcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    %scatter(X,get_ax_val(.35,ylim),Msz,Rcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    X = get_ax_val(X0+(plt_idx*move_pt),xlim);
    
   % %regular w/ colors
   % scatter(X,get_ax_val(Y0-.05,ylim),Msz,Bcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
   % scatter(X,get_ax_val(Y0-.15,ylim),Msz,Rcol(col_idx,:),'filled','MarkerFaceAlpha',alph)
    
    %for symbols
    scatter(X,get_ax_val(Y0-.05,ylim),Msz,Bcol(col_idx,:),'Marker',markers{plt_idx})
    scatter(X,get_ax_val(Y0-.15,ylim),Msz,Rcol(col_idx,:),'Marker',markers{plt_idx})
end

set(gca,'FontSize',18)
switch print_anything
    case 'yes'
        print(fullfile(figdir,[fig_fn '-scatter']),'-djpeg')
        savefig(fullfile(figdir,[fig_fn '-scatter']))
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


