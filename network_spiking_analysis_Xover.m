clear
clc
hold off; close all
format compact

num_workers = 24;

sim_name = 'network_spiking';
name4plots = 'Xover'; %added to basic fig directory, distingish different timecourse plots etc
home_dir = '~/Desktop/ksander/rotation/project';
addpath(home_dir)
addpath(fullfile(home_dir,'helper_functions'))
fig_dir = fullfile(home_dir,'Results','network_spiking_figures');
fig_dir = fullfile(fig_dir,name4plots,'cellwise');
resdir = fullfile(home_dir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);

Xgroups = {'I-stay','E-stay'};
Tcourse = 'presw250to5'; %'preswitch' | 'all' |presw250to5
treat_data = 'base0'; % zscore | base0 | minmax
Nstraps = 25e3;
print_anything = 'yes';

c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('intersections.m')})

switch treat_data
    case 'zscore'
        fig_dir = fullfile(fig_dir,'zscore');
    case 'base0'
        fig_dir = fullfile(fig_dir,'base0');
    case 'minmax'
        fig_dir = fullfile(fig_dir,'minmax');
end
if ~isdir(fig_dir),mkdir(fig_dir);end


switch Tcourse
    case 'preswitch'
        fig_fn = 'preswitch_timecourse';
        preswitch_win = 105e-3; %preswitch window to consider (T0-X)
        postswitch_win = -5e-3; %postwitch window to consider (T+X)
    case 'all'
        fig_fn = 'switching_timecourse';
        preswitch_win = 250e-3; %preswitch window to consider (T0-X)
        postswitch_win = 150e-3; %postwitch window to consider T+X)
    case 'presw250to5'
        fig_fn = 'presw250to5';
        preswitch_win = 250e-3; %preswitch duration to plot (T0-X)
        postswitch_win = -5e-3; %postwitch duration to plot (T+
end
fig_fn = sprintf('%s_%sx%s',fig_fn,Xgroups{:});
fig_fn = strrep(fig_fn,'-','');


recorded_switchtime =  250e-3; %actual switchtime in recorded switch
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%result summaries
fontsz = 30;
lnsz = 3; %spikerate plots
orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];
%IMPORTANT: legend_labels is used for ordering stuff in this analysis tooo
legend_labels = {'E-stay','E-switch','I-stay','I-switch'};


%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);
%match this to the legend labels
celltype.Estay = celltype.excit & celltype.pool_stay;
celltype.Eswitch = celltype.excit & celltype.pool_switch;
celltype.Istay = celltype.inhib & celltype.pool_stay;
celltype.Iswitch = celltype.inhib & celltype.pool_switch;

%replace with like... "network types"
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

%get results & sort into network type
num_files = numel(output_fns);
num_types = numel(network_pair_info);
%reformat the parameter sets for easier searching/indexing
Psets = num2cell(cat(1,network_pair_info{:}),2);
%file_data = cell(num_files,2);
result_data = num2cell(zeros(size(Psets)));
switch_counts = zeros(size(result_data));

% for idx = 1:num_files
%     if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
%     curr_file = load(output_fns{idx});
%     %find parameter set
%     p = curr_file.options;
%     p = {p.ItoE,p.EtoI,unique(p.trial_stimuli),p.stim_targs};
%     p = cellfun(@(x) isequal(x,p),Psets); %index
%     %file data
%     fdata = curr_file.sim_results{2};
%     switch_counts(p) = switch_counts(p) + size(fdata,3); %keep track of how many timecourses
%     fdata = sum(fdata,3); %sum over the switches from this file
%     %add to the rest of them
%     result_data{p} = result_data{p} + fdata;
% end
% 
% %very quickly check for guys that don't have many observations
% for idx = 1:numel(switch_counts)
%     if switch_counts(idx) < 10000
%         fprintf('\nnetwork:\n--- [ItoE: %.2f] [EtoI: %.2f] [stim: %.2f hz] [%s]\nhas < 10k states (%i)\n',...
%             Psets{idx}{:},switch_counts(idx))
%     end
% end
% %save('checkpoint','result_data','switch_counts')
% %now divide the summed spike timecourses
% switch_counts = num2cell(switch_counts);
% result_data = cellfun(@(x,y) x./y,result_data,switch_counts,'UniformOutput',false);
% result_data = cellfun(@(x) sim_spikerate(x,timestep),result_data,'UniformOutput',false);


load('checkpoint.mat')


%only consider -Xms to +Xms
recorded_switchtime = recorded_switchtime/timestep; %actual switchtime in recorded switch
postswitch_win = postswitch_win/timestep;
preswitch_win = preswitch_win/timestep;
%record_duration = size(cat(1,result_data{:}),2); %get the duration of recorded timecourses
data_window = 1 + recorded_switchtime - preswitch_win:recorded_switchtime + postswitch_win;
onset_switch = 1 + recorded_switchtime - min(data_window); %adjusted to the new plotting window
%cut down the data matrix to this window
result_data = cellfun(@(x) x(:,data_window,:),result_data,'UniformOutput',false);

switch treat_data
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

%nah, this time overlap the plots--- so it's 5x1 subplot
title_lab = sprintf('%s x %s',Xgroups{:});

%Xax = ((data_window-onset_switch).*timestep)./1e-3;
orient tall

for idx = 1:num_types
    disp(fprintf('working on network #%i / 5 ...\n',idx))
    h(idx) = subplot(5,1,idx);
    curr_net_info = network_pair_info{idx};
    hold on
    Xtime = NaN(1,2);
    Xnull = NaN(Nstraps,2);
    for j = 1:2
        
        %get the right color
        if strcmpi(curr_net_info{j,end},'Eswitch')
            Nspeed = 'slow';
            lcol = orange;
        elseif strcmpi(curr_net_info{j,end},'Estay')
            Nspeed = 'fast';
            lcol = matblue;
        end
        
        %find the right results for network set-up
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),Psets,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = result_data{curr_data};
        Xt = crossover_analysis(curr_data,celltype,Xgroups);
        %returns x,y intersection point
        %actual time in ms
        Xt = Xt(:,1);
        Xt = ((Xt-onset_switch).*timestep) ./ 1e-3;
        Xtime(j) = Xt;
        
        %bootstrap crossing distribution
        X0 = crossover_analysis(curr_data,celltype,Xgroups,Nstraps,'parallel');
        X0 = X0(:,1);
        X0 = ((X0-onset_switch).*timestep) ./ 1e-3;
        Xnull(:,j) = X0;
        d(j) = histogram(X0,'FaceColor',lcol,'EdgeColor',lcol,'BinWidth',1); %1ms bins
        
        
    end
    
    X_yvals = get(gca,'YLim');
    X_yvals = linspace(X_yvals(1),X_yvals(2));
    for pl_idx = 1:2
        %plot the true means nicely scaled to Y axis
        plot(repmat(Xtime(pl_idx),1,numel(X_yvals)),X_yvals,'LineWidth',lnsz,'Color',d(pl_idx).FaceColor)
        scatter(Xtime(pl_idx),X_yvals(end),50,'k','*')
    end
    
    %set(gca,'XLim',[Xax(1),Xax(end)]);
    %Xticks = num2cell(get(gca,'Xtick'));
    %Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
    %set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    leg_info = curr_net_info(:,end)';
    leg_info = strrep(leg_info,'Eswitch','slow');
    leg_info = strrep(leg_info,'Estay','fast');
    legend(leg_info,'location','best')
    hold off
    if idx == num_types
        xlabel('Leave decision (ms)')
    end
    if idx == 1
        title(title_lab,'Fontsize',14)
    end
    ylabel(sprintf('net #%i %s',idx))
    axis tight
    
end
linkaxes(h,'x')
%make the tick labels nice 
for idx = 1:num_types
    Xticks = num2cell(get(h(idx),'Xtick'));
    Xlabs = cellfun(@(x) sprintf('%+i',round(x)),Xticks,'UniformOutput', false); %this is for normal stuff
    set(h(idx), 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
end

switch print_anything
    case 'yes'
        print(fullfile(fig_dir,fig_fn),'-djpeg','-r300')
end

delete(gcp('nocreate'))



function cell_raster = sim_spikerate(cell_raster,timestep)

num_binsamps = 2e-3/timestep; %num samples in 2ms
raster_sz = size(cell_raster);
bin_magic = [numel(cell_raster(:,1)), num_binsamps, numel(cell_raster(1,:))/num_binsamps]; %set up for a magic trick
cell_raster = reshape(cell_raster,bin_magic);
cell_raster = squeeze(sum(cell_raster,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_raster = repmat(cell_raster,[num_binsamps 1 1]);
cell_raster = reshape(cell_raster,raster_sz); %put the rabit back in the hat
end




