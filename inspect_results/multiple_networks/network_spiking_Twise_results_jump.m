clear
clc
format compact
hold off;close all

%note-- 4/22/19: check when, for a specific cell pool, there's a large jump
%in activity. For instance, when does the spikereate increase by 5x


opt = struct();
opt.multiple_stimuli = 'no'; %'yes'|'no';
opt.params2match = {'conn','stim'}; %specify how results are matched to network types (at most {'conn','stim'})
opt.print_anything = 'yes'; %'yes' | 'no';
opt.Tcourse = 'all'; %'preswitch' | 'all' | 'presw250to5' | 'presw150to25'
opt.treat_data = 'none'; %'base0'; % zscore | base0 | minmax | 'none'
opt.Ypool = 'I-stay'; %Y data for cross correlation
opt.jump_factor = 5; %5x increase
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)

%specify simulation
%---sim setup-----------------

Snames = {'nets_slowD','nets_slowD_baseline'}; %runs outa memory with basline data
figdir = {'figures_nets_slowD','figures_nets_slowD'};

basedir = '~/Desktop/ksander/rotation/project/';
addpath(fullfile(basedir,'helper_functions'))

%loop over these
%timewins = {'preswitch', 'all', 'presw250to5', 'presw150to25'};
timewins = {'presw250to5'};
treatments = {'none'}; %{'none', 'base0', 'minmax'};
Yjump = 2:2:10;

for Sidx = 1:numel(Snames)
    
    %loop through & do everything for these results
    for Tidx = 1:numel(timewins)
        
        opt.Tcourse = timewins{Tidx};
        
        for Didx = 1:numel(treatments)
            
            opt.treat_data = treatments{Didx};
            
            for Jidx = 1:numel(Yjump)
                
                opt.jump_factor = Yjump(Jidx);
                
                %make the figures
                netspiking_figure(basedir,Snames{Sidx},figdir{Sidx},opt)
            end
        end
    end
end





function netspiking_figure(home_dir,sim_name,figdir,opt)
hold off;close all

mainfig_dir = fullfile(home_dir,'Results',figdir,'Twise_jump');
resdir = fullfile(home_dir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
params2match = opt.params2match;
fig_dir = fullfile(mainfig_dir,sprintf('jump_factor_%i',opt.jump_factor)); %strrep(opt.treat_data,'none','no_treatment'));
data_fn = sprintf('data_%s_%s_%s.mat',sim_name,opt.Ypool,opt.Tcourse);
if exist(fullfile(fig_dir,data_fn)) > 0,load_summary = true;else,load_summary = false;end
if contains(sim_name,'baseline'),BLdata = true;else,BLdata = false;end

recorded_switchtime =  250e-3; %actual switchtime in recorded switch

CIalpha = .05;
rate_binsz = 20e-3; %binsize for spikerate calculations

switch opt.pulse_stim
    case 'off'
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        error('get this from  options dude, was previously striped from sim name')
end

switch opt.Tcourse
    case 'preswitch'
        fig_fn = 'preswitch_timecourse';
        preswitch_plottime = 105e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = -5e-3; %postwitch duration to plot (T+X)
    case 'all'
        fig_fn = 'switching_timecourse';
        preswitch_plottime = 250e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = 150e-3; %postwitch duration to plot (T+X)
    case 'presw250to5'
        fig_fn = 'presw250to5';
        preswitch_plottime = 250e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = -5e-3; %postwitch duration to plot (T+X)
    case 'presw150to25'
        fig_fn = 'presw150to25';
        preswitch_plottime = 150e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = 25e-3; %postwitch duration to plot (T+X)
end

fig_fn = sprintf('%s_%s_%s',sim_name,opt.Ypool,fig_fn);
if ~isdir(fig_dir),mkdir(fig_dir);end

Yax_labs = 'p(x)';

%result summaries
fontsz = 30;
lnsz = 3; %spikerate plots
orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];
legend_labels = {'E-stay','E-switch','I-stay','I-switch'};
figure;hold on;lncols = []; %get color order
for idx = 1:numel(legend_labels)
    lncols(idx) = plot(NaN,'Linewidth',lnsz);
end
lncols = get(lncols,'Color');hold off;close all
cellcat = @(x,y) cat(y,x{:}); %faster cell2mat...

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);

%make a pool average function, match to legend labels ordering
pool_inds = cell(size(legend_labels))';
pool_inds{1} = celltype.excit & celltype.pool_stay; %E-stay
pool_inds{2} = celltype.excit & celltype.pool_switch; %E-switch
pool_inds{3} = celltype.inhib & celltype.pool_stay; %I-stay
pool_inds{4} = celltype.inhib & celltype.pool_switch; %I-switch
%function for mean pool timepoint. Takes 2D rasters, gives cell
pool_means = @(d) cellfun(@(x) mean(d(x,:),1),pool_inds,'UniformOutput',false);

%indicies for cross-correlation
Ydata_inds = strcmpi(legend_labels',opt.Ypool);
Xdata_inds = ~Ydata_inds;

%get general options file from the first file
gen_options = load(output_fns{1});
gen_options = gen_options.options;
timestep = gen_options.timestep;

switch opt.multiple_stimuli
    case 'yes'
        param_varnams = {'ItoE','EtoI','stim_A','stim_B','targ_cells'};
        error('not configured yet'); %look at line below, figure out what you gotta do
        %IDvars = param_varnams(~ismember(param_varnams,'stim_B')); %stim B not particular to network type
    case 'no'
        param_varnams = {'ItoE','EtoI','stim','targ_cells'};
end
%for indexing the result paramters
IDvars = [];
if sum(strcmp('conn',params2match)) > 0,IDvars = {'ItoE','EtoI'};end
if sum(strcmp('stim',params2match)) > 0,IDvars = [IDvars,param_varnams(startsWith(param_varnams,'stim'))];end


%info on the specific network parameters in this simulation
num_net_types = 10;
num_pairs = 5;
pair_inds = num2cell(reshape(1:num_net_types,[],num_pairs)); %just gives a cell array for pair indicies
network_pair_info = cell(num_pairs,1);
Psets = [];
for idx = 1:num_pairs
    curr_params = cellfun(@(x) get_network_params(x,gen_options),pair_inds(:,idx),'UniformOutput',false);
    switch opt.multiple_stimuli
        case 'yes'
            curr_params = cellfun(@(x)...
                {x.ItoE, x.EtoI,x.trial_stimuli,x.stim_targs},...
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
    if BLdata %change net paramters to reflect baseline sim
        T{:,param_varnams(startsWith(param_varnams,'stim'))} = 0;
        T.targ_cells(:) = {'baseline'};
    end
    network_pair_info{idx} = T;
    Psets = [Psets;table2cell(T)];
end
Psets = num2cell(Psets,2);


%only plot -Xms to +Xms (need round() for roundoff errors...)
recorded_switchtime = round(recorded_switchtime/timestep,4); %actual switchtime in recorded switch
postswitch_plottime = round(postswitch_plottime/timestep,4);
preswitch_plottime = round(preswitch_plottime/timestep,4);
%record_duration = size(cat(1,result_data{:}),2); %get the duration of recorded timecourses
plotting_window = 1 + recorded_switchtime - preswitch_plottime:recorded_switchtime + postswitch_plottime;
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window
%onset_switch = onset_switch - round(gen_options.state_test_time ./ timestep,4); %adjust for the threshold testing time
%cut down the data matrix to this window
PWwin_start = plotting_window(1); PWwin_stop = plotting_window(end);
cutdown_data = @(x) x(:,PWwin_start:PWwin_stop,:);

switch opt.treat_data
    case 'base0'
        data_treatment = @(D) cellfun(@(x) bsxfun(@minus, x,  min(x,[],2)),D,'UniformOutput',false);
    case 'minmax'
        minmax_norm = @(x) bsxfun(@rdivide, bsxfun(@minus,x,min(x,[],2)) ,max(x,[],2) - min(x,[],2));
        data_treatment = @(D) cellfun(minmax_norm,D,'UniformOutput',false);
    case 'none'%parpool gets mad if this isn't defined
        data_treatment = @(D) cellfun(@(x) x ,D,'UniformOutput',false);
end

num_files = numel(output_fns);
file_data = cell(num_files,2);

if ~load_summary
    
    num_workers = 24;
    c = parcluster('local');
    c.NumWorkers = num_workers;
    parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_stay_durations'),which('xcorr')})
    special_progress_tracker = fullfile(home_dir,'SPT.txt');
    if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start
    
    parfor idx = 1:num_files
        curr_file = load(output_fns{idx});
        %verify recorded spiking results are valid...
        all_events = curr_file.sim_results{1};
        [~,valid_events] = find_stay_durations(all_events,curr_file.options,'verify');
        valid_events = cat(1,valid_events{:,1});
        fevents = curr_file.sim_results{3}(:,1); %time indicies for the recorded spiking events
        fevents = cat(1,fevents{:}) .* curr_file.options.timestep; %convert to time for comparison
        valid_events = ismember(fevents,valid_events);
        %file data
        fdata = curr_file.sim_results{2};
        fdata = fdata(:,:,valid_events); %ensure records are valid
        fdata = squeeze(num2cell(fdata,1:2)); %unclear why num2cell(x,3) didn't do this..
        fdata = cellfun(pool_means,fdata,'UniformOutput',false); %mean pool spikes per timepoint
        fdata = cellfun(@(x) cellcat(x,1),fdata,'UniformOutput',false);
        fdata = cellfun(@(x) sim_spikerate(x,timestep,rate_binsz),fdata,'UniformOutput',false); %spikerate per bin
        %apply the data treatment
        switch opt.treat_data
            case 'none'
            otherwise
                fdata = data_treatment(fdata);
        end
        fdata = cellfun(@(x) x(Ydata_inds,:),fdata,'UniformOutput',false); %only the pool of interest
        fdata = cellfun(cutdown_data,fdata,'UniformOutput',false); %to specified analysis window
        
        for Ti = 1:numel(fdata)
            curr_trace = fdata{Ti};
            curr_trace(curr_trace <1) = 1; %set 1hz min
            curr_trace = curr_trace(2:end) ./ curr_trace(1:end-1);
            curr_trace = [1,curr_trace];
            curr_trace = find(curr_trace >= opt.jump_factor,1,'first');
            curr_trace = (curr_trace - onset_switch)  * timestep;
            fdata{Ti} = curr_trace;
        end
        fdata = fdata(~cellfun(@isempty,fdata));
        fdata = cellcat(fdata,1);
        %store results & parameters
        file_data(idx,:) = {fdata,curr_file.options};
        progress = worker_progress_tracker(special_progress_tracker);
        if mod(progress,floor(num_files * .05)) == 0 %at 5 percent
            progress = (progress / num_files) * 100;
            fprintf('%s --- %.1f percent complete\n',datestr(now,31),progress);
        end
    end
    
    delete(gcp('nocreate'))
    delete(special_progress_tracker)
    
    fprintf('\nsaving data...\n')
    save(fullfile(fig_dir,data_fn),'file_data')
    
elseif load_summary
    fprintf('\nloading saved summary data...\n')
    file_data = load(fullfile(fig_dir,data_fn));
    file_data = file_data.file_data;
end



%collapse jobs
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'};
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
    result_data{idx,1} = cellcat(file_data(curr_file,1),1);
    fprintf('------n states = %i\n',size(result_data{idx,1},3))
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
    
    %delete this so you don't run outa memory
    file_data(curr_file,1) = cell(size(file_data(curr_file,1)));
end

clear file_data

%this is stupid & obviously a hold-over from something I didn't implement well in the first place...
net_type.targ_cells = cellfun(@(x) stimtarg_labels{x},...
    num2cell(net_type.targ_cells),'UniformOutput',false);


ln.baseline.Color = [103 115 122] ./ 255;
plt_idx = 0;
for idx = 1:num_pairs
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        
        %get the right color
        Nspeed = curr_net_info.Row{j};
        Ntargs = curr_net_info.targ_cells{j};
        
        %find the right results for network set-up
        curr_inds = cellfun(@(x) isequal(x,table2cell(curr_net_info(j,:))),Psets,'UniformOutput',false);
        curr_inds = cat(1,curr_inds{:});
        curr_data = result_data{curr_inds};
        curr_data = curr_data ./ 1e-3;
        histogram(curr_data,'Normalization','probability');
        %legend({sprintf('\\mu = %.0f\nMed = %.0f',mean(curr_data),median(curr_data))},...
        %    'Location','best','Box','off')
        text(1,1,sprintf('\\mu = %.0f \nMed = %.0f ',mean(curr_data),median(curr_data)),...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
            
        hold off
        if plt_idx == 9 || plt_idx == 10
            xlabel(sprintf('T - state transition (ms)\n%s rate increases by %ix',opt.Ypool,opt.jump_factor))
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',Nspeed),'Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('net #%i\n%s',idx,Yax_labs))
        end
    end
end
orient tall

axis tight

linkaxes(h,'xy')

% fake_ax = axes('Position',[-1,-1,0,0],'Visible','off');
% Xcolors =lncols(Xdata_inds);
% hold on
% for idx = 1:numel(legend_labels(Xdata_inds))
%     lns(idx) = plot(NaN,'Linewidth',lnsz,'Color',Xcolors{idx},'Parent',fake_ax);
% end
% hold off
% lp = legend(lns,legend_labels(Xdata_inds),'FontWeight','b','Fontsize',14,...
%     'Location','northoutside','Box','off','Orientation','horizontal');
% lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];

switch opt.print_anything
    case 'yes'
        print(fullfile(fig_dir,fig_fn),'-djpeg','-r300')
end


close all;hold off

end