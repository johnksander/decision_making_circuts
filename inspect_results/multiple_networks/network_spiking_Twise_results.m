clear
clc
format compact
hold off;close all

%note-- 4/16/19: this is the code for analyzing trial-wise simulation spikerates

%to do:

opt = struct();
opt.multiple_stimuli = 'no'; %'yes'|'no';
opt.params2match = {'conn','stim'}; %specify how results are matched to network types (at most {'conn','stim'})
opt.print_anything = 'yes'; %'yes' | 'no';
opt.Tcourse = 'all'; %'preswitch' | 'all' | 'presw250to5' | 'presw150to25'
opt.treat_data = 'none'; %'base0'; % zscore | base0 | minmax | 'none'
opt.zoomed_fig = 'no'; %'yes'|'no'; %ignore I-stay spiking for Y limits
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)

%specify simulation
%---sim setup-----------------

Snames = {'nets_fastD','nets_fastD_baseline'};
figdir = {'figures_nets_fastD'};

basedir = '~/Desktop/work/ACClab/rotation/project/'; %'~/Desktop/ksander/rotation/project/';
addpath(fullfile(basedir,'helper_functions'))

%loop over these
timewins = {'preswitch', 'all', 'presw250to5', 'presw150to25'};
treatments = {'none','zscore', 'base0', 'minmax'};
zooming = {'yes', 'no'};

for Sidx = 1:numel(Snames)
    
    %loop through & do everything for these results
    for Tidx = 1:numel(timewins)
        
        opt.Tcourse = timewins{Tidx};
        
        for Didx = 1:numel(treatments)
            
            opt.treat_data = treatments{Didx};
            
            for Zidx  = 1:numel(zooming)
                
                opt.zoomed_fig = zooming{Zidx};
                
                %make the digures
                netspiking_figure(basedir,Snames{Sidx},figdir{Sidx},opt)
            end
        end
    end
end





function netspiking_figure(home_dir,sim_name,figdir,opt)
hold off;close all

mainfig_dir = fullfile(home_dir,'Results',figdir,'spikeplots_Twise');
resdir = fullfile(home_dir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
params2match = opt.params2match;
data_fn = sprintf('summary_Twise_data_%s.mat',sim_name);
if exist(fullfile(mainfig_dir,data_fn)) > 0,load_summary = true;else,load_summary = false;end
if contains(sim_name,'baseline'),BLdata = true;else,BLdata = false;end


switch opt.pulse_stim
    case 'off'
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        error('get this from  options dude, was previously striped from sim name')
end

recorded_switchtime =  250e-3; %actual switchtime in recorded switch
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

fig_fn = sprintf('%s_%s',sim_name,fig_fn);

switch opt.treat_data
    case 'zscore'
        fig_dir = fullfile(mainfig_dir,'zscore');
        fig_fn = sprintf('%s_%s',fig_fn,'zscore');
        Yax_labs = 'Z spiking';
    case 'base0'
        fig_dir = fullfile(mainfig_dir,'base0');
        fig_fn = sprintf('%s_%s',fig_fn,'base0');
        Yax_labs = 'rate - min (Hz)';
    case 'minmax'
        fig_dir = fullfile(mainfig_dir,'minmax');
        fig_fn = sprintf('%s_%s',fig_fn,'minmax');
        Yax_labs = 'spiking (0-1)';
    case 'none'
        fig_dir = fullfile(mainfig_dir,'no_treatment');
        Yax_labs = 'spiking (Hz)';
end

switch opt.zoomed_fig
    case 'yes'
        fig_dir = fullfile(fig_dir,'zoomed');
        fig_fn = [fig_fn '_zoomed'];
    case 'no'
        axis tight
end

if ~isdir(fig_dir),mkdir(fig_dir);end

%result summaries
fontsz = 30;
lnsz = 3; %spikerate plots
orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];
legend_labels = {'E-stay','E-switch','I-stay','I-switch'};


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

pool_means = @(x) 
        
    
        Estay = mean(curr_data(celltype.excit & celltype.pool_stay,:),1);
        Eswitch = mean(curr_data(celltype.excit & celltype.pool_switch,:),1);
        Istay = mean(curr_data(celltype.inhib & celltype.pool_stay,:),1);
        Iswitch = mean(curr_data(celltype.inhib & celltype.pool_switch,:),1);
        


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


if ~load_summary
    %get results & sort into network type
    num_files = numel(output_fns);
    result_data = num2cell(zeros(size(Psets)));
    switch_counts = zeros(size(result_data));
    
    for idx = 1:num_files
        if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
        curr_file = load(output_fns{idx});
        %find parameter set
        p = curr_file.options;
        switch opt.multiple_stimuli
            case 'yes'
                error('not configed')
            otherwise
                p = {p.ItoE,p.EtoI,unique(p.trial_stimuli),p.stim_targs};
        end
        p = cellfun(@(x) isequal(x,p),Psets); %index
        
        %verify recorded spiking results are valid... after-the-fact... this
        %will be fixed in the recording paramters for the next simulation job
        all_events = curr_file.sim_results{1};
        [~,valid_events] = find_stay_durations(all_events,curr_file.options,'verify');
        valid_events = cat(1,valid_events{:,1});
        fevents = curr_file.sim_results{3}(:,1); %time indicies for the recorded spiking events
        fevents = cat(1,fevents{:}) .* curr_file.options.timestep; %convert to time for comparison
        valid_events = ismember(fevents,valid_events);
        %now back to more regular stuff
        keyboard
        %file data
        fdata = curr_file.sim_results{2};
        fdata = fdata(:,:,valid_events); %ensure records are valid
        
        
        result_data = cellfun(@(x,y) x./y,result_data,switch_counts,'UniformOutput',false);
        result_data = cellfun(@(x) sim_spikerate(x,timestep),result_data,'UniformOutput',false);
        
        Estay = mean(curr_data(celltype.excit & celltype.pool_stay,:),1);
        Eswitch = mean(curr_data(celltype.excit & celltype.pool_switch,:),1);
        Istay = mean(curr_data(celltype.inhib & celltype.pool_stay,:),1);
        Iswitch = mean(curr_data(celltype.inhib & celltype.pool_switch,:),1);
        
        
        switch_counts(p) = switch_counts(p) + size(fdata,3); %keep track of how many timecourses
        fdata = sum(fdata,3); %sum over the switches from this file
        %add to the rest of them
        result_data{p} = result_data{p} + fdata;
    end

    %now divide the summed spike timecourses
    switch_counts = num2cell(switch_counts);
    result_data = cellfun(@(x,y) x./y,result_data,switch_counts,'UniformOutput',false);
    result_data = cellfun(@(x) sim_spikerate(x,timestep),result_data,'UniformOutput',false);
    
    fprintf('\nsaving data...\n')
    save(fullfile(mainfig_dir,data_fn),'result_data','switch_counts')
elseif load_summary
    fprintf('\nloading saved summary data...\n')
    Fdata = load(fullfile(mainfig_dir,data_fn));
    result_data = Fdata.result_data;
    switch_counts = Fdata.switch_counts;
end

%very quickly check for guys that don't have many observations
for idx = 1:numel(switch_counts)
    if switch_counts{idx} < 10000
        fprintf('\nnetwork:\n--- [ItoE: %.2f] [EtoI: %.2f] [stim: %.2f hz] [%s]\nhas < 10k states (%i)\n',...
            Psets{idx}{:},switch_counts{idx})
    end
end


%only plot -Xms to +Xms (need round() for roundoff errors...)
recorded_switchtime = round(recorded_switchtime/timestep,4); %actual switchtime in recorded switch
postswitch_plottime = round(postswitch_plottime/timestep,4);
preswitch_plottime = round(preswitch_plottime/timestep,4);
%record_duration = size(cat(1,result_data{:}),2); %get the duration of recorded timecourses
plotting_window = 1 + recorded_switchtime - preswitch_plottime:recorded_switchtime + postswitch_plottime;
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window
%cut down the data matrix to this window
result_data = cellfun(@(x) x(:,plotting_window,:),result_data,'UniformOutput',false);

switch opt.treat_data
    case 'zscore'
        result_data = cellfun(@(x) zscore(x,[],2),result_data,'UniformOutput',false);
    case 'base0'
        result_data = cellfun(@(x) bsxfun(@minus, x,  min(x,[],2)),result_data,'UniformOutput',false);
    case 'minmax'
        Xmin = cellfun(@(x) min(x,[],2),result_data,'UniformOutput',false);
        Xmax = cellfun(@(x) max(x,[],2),result_data,'UniformOutput',false);
        result_data = cellfun(@(x,y) bsxfun(@minus,x,y),result_data,Xmin,'UniformOutput',false);
        result_data = cellfun(@(x,y,z) bsxfun(@rdivide,x,z-y),result_data,Xmin,Xmax,'UniformOutput',false);
end

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
        curr_data = cellfun(@(x) isequal(x,table2cell(curr_net_info(j,:))),Psets,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = result_data{curr_data};
        
        %plot by celltype
        %normal people indexing that makes sense, then take the mean
        Estay = mean(curr_data(celltype.excit & celltype.pool_stay,:),1);
        Eswitch = mean(curr_data(celltype.excit & celltype.pool_switch,:),1);
        Istay = mean(curr_data(celltype.inhib & celltype.pool_stay,:),1);
        Iswitch = mean(curr_data(celltype.inhib & celltype.pool_switch,:),1);
        
        
        %plot the aggregated timecourses
        hold on
        ln.Estay = plot(Estay,'Linewidth',lnsz);
        ln.Eswitch = plot(Eswitch,'Linewidth',lnsz);
        plot(Istay,'Linewidth',lnsz)
        plot(Iswitch,'Linewidth',lnsz)
        
        leg_col = ln.(Ntargs).Color;
        l = plot(NaN,'Color',leg_col,'LineWidth',lnsz); %for legend
        
        hold off
        Xlim_max = numel(curr_data(1,:));
        set(gca,'XLim',[0 Xlim_max]);
        set(gca,'Xtick',linspace(0,Xlim_max,5));
        Xticks = num2cell(get(gca,'Xtick'));
        if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
        Xlabs = cellfun(@(x) sprintf('%+i',round(((x-onset_switch)*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        switch opt.pulse_stim
            case 'off'
                leg_labels = sprintf('%.0f Hz %s',curr_net_info.stim(j),curr_net_info.targ_cells{j});
                l = legend(l,leg_labels,'location','best','Box','off');l.Box = 'off'; %Why do I need boxoff 2x??
            otherwise
                legend(l,sprintf('%.0fHz %s',curr_net_info{j,3},Pdur),'location','best','box','off')
        end
        warning('on','MATLAB:legend:IgnoringExtraEntries')
        
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('Leave decision (ms)')
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

switch opt.zoomed_fig
    case 'yes'
        %zoom in better
        Yl = arrayfun(@(x) x.Children,h,'UniformOutput',false);
        %skip legend & I-stay
        keep_inds = find(~ismember(fliplr(legend_labels),'I-stay')) + 1;
        Yl = cellfun(@(x) x(keep_inds),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) cat(1,x(:).YData),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) max(x(:)),Yl);
        arrayfun(@(x,y) ylim(x,[0,y]),h,Yl,'UniformOutput',false);
    case 'no'
        axis tight
end

linkaxes(h,'x')

fake_ax = axes('Position',[-1,-1,0,0],'Visible','off');
hold on
for idx = 1:numel(legend_labels)
    lns(idx) = plot(NaN,'Linewidth',lnsz,'Parent',fake_ax);
end
hold off
lp = legend(lns,legend_labels,'FontWeight','b','Fontsize',14,...
    'Location','northoutside','Box','off','Orientation','horizontal');
lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];


switch opt.print_anything
    case 'yes'
        print(fullfile(fig_dir,fig_fn),'-djpeg','-r300')
end


close all;hold off
end