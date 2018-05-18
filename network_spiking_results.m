clear
clc
hold off; close all
format compact

%this comes from writeup_inspect_results_combfig_20170725.m
%changing the figures a bit for brownbag
%NOTE: this script requires 2017a. Local function below.
%home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
sim_name = 'network_spiking';
name4plots = 'combfig'; %added to basic fig directory, distingish different timecourse plots etc
home_dir = '~/Desktop/ksander/rotation/project';
addpath(home_dir)
addpath(fullfile(home_dir,'helper_functions'))
fig_dir = fullfile(home_dir,'Results','network_spiking_figures');
fig_dir = fullfile(fig_dir,name4plots);
resdir = fullfile(home_dir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
if ~isdir(fig_dir),mkdir(fig_dir);end
fig_name = 'combfig';
Tcourse = 'preswitch'; %'preswitch' | 'all'
print_anything = 'yes';
zoomed_fig = 'no'; %ignore I-stay spiking for Y limits

switch Tcourse
    case 'preswitch'
        fig_fn = 'preswitch_timecourse';
        preswitch_plottime = 105e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = -5e-3; %postwitch duration to plot (T+X)
    case 'all'
        fig_fn = 'switching_timecourse';
        preswitch_plottime = 250e-3; %preswitch duration to plot (T0-X)
        postswitch_plottime = 150e-3; %postwitch duration to plot (T+X)
end

recorded_switchtime =  250e-3; %actual switchtime in recorded switch
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

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
% %very quickly check for guys that don't have many observations
% for idx = 1:numel(switch_counts)
%     if switch_counts(idx) < 10000
%         fprintf('\nnetwork:\n--- [ItoE: %.2f] [EtoI: %.2f] [stim: %.2f hz] [%s]\nhas < 10k states (%i)\n',...
%             Psets{idx}{:},switch_counts(idx))
%     end
% end
%save('checkpoint','result_data','switch_counts')
load('checkpoint.mat')

%now divide the summed spike timecourses 
switch_counts = num2cell(switch_counts);
result_data = cellfun(@(x,y) x./y,result_data,switch_counts,'UniformOutput',false);
result_data = cellfun(@(x) sim_spikerate(x,timestep),result_data,'UniformOutput',false);

figIDs = 'spikes';
legloc = 'west';
Yax_labs = 'spike rate (Hz)';

%only plot -Xms to +Xms
recorded_switchtime = recorded_switchtime/timestep; %actual switchtime in recorded switch
postswitch_plottime = postswitch_plottime/timestep;
preswitch_plottime = preswitch_plottime/timestep;
record_duration = size(cat(1,result_data{:}),2); %get the duration of recorded timecourses
plotting_window = 1 + recorded_switchtime - preswitch_plottime:recorded_switchtime + postswitch_plottime;
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window
%cut down the data matrix to this window
result_data = cellfun(@(x) x(:,plotting_window,:),result_data,'UniformOutput',false);


plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        
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
                
        %plot by celltype
        
        %normal people indexing that makes sense, then take the mean 
        Estay = mean(curr_data(celltype.excit & celltype.pool_stay,:),1);
        Eswitch = mean(curr_data(celltype.excit & celltype.pool_switch,:),1);
        Istay = mean(curr_data(celltype.inhib & celltype.pool_stay,:),1);
        Iswitch = mean(curr_data(celltype.inhib & celltype.pool_switch,:),1);
        
          
        %plot the aggregated timecourses
        hold on
        plot(Estay,'Linewidth',lnsz)
        plot(Eswitch,'Linewidth',lnsz)
        plot(Istay,'Linewidth',lnsz)
        plot(Iswitch,'Linewidth',lnsz)
        l = plot(NaN,'Color',lcol,'LineWidth',lnsz); %for legend
        
        hold off
        Xlim_max = numel(curr_data(1,:));
        set(gca,'XLim',[0 Xlim_max]);
        Xticks = num2cell(get(gca,'Xtick'));
        %Xticks = Xticks(2:2:end-1);
        if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
        Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        legend(l,sprintf('%.0fHz',curr_net_info{j,3}),'location','best')
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('Leave decision (ms)')
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',Nspeed),'Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('net #%i spiking (hz)',idx))
        end
          
    end
end
orient tall

switch zoomed_fig
    case 'yes'
        %zoom in better
        Yl = arrayfun(@(x) x.Children,h,'UniformOutput',false);
        %skip legend & I-stay
        Yl = cellfun(@(x) x([2,4,5]),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) cat(1,x(:).YData),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) max(x(:)),Yl);
        arrayfun(@(x,y) ylim(x,[0,y]),h,Yl,'UniformOutput',false);
        fig_fn = [fig_fn '_zoomed'];
    case 'no'
        axis tight
end

linkaxes(h,'x')
switch print_anything
    case 'yes'
        print(fullfile(fig_dir,fig_fn),'-djpeg','-r300')
end

close all;hold off
%print a legend seperately

%plot the aggregated timecourses
hold on
for idx = 1:numel(legend_labels)
    tr_col(idx) = plot(NaN,'Linewidth',lnsz);
end
hold off
legend(legend_labels)
switch print_anything
    case 'yes'
        print(fullfile(fig_dir,sprintf('%s_legend',fig_fn)),'-djpeg','-r300')
end

%save the line colors for below 
tr_col = arrayfun(@(x) x.Color,tr_col,'UniformOutput',false);
        
%do the like... crossing time plot 
close all;hold off


plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        
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
                
        %plot by celltype
        
        %normal people indexing that makes sense, then take the mean 
        Estay = mean(curr_data(celltype.excit & celltype.pool_stay,:),1);
        Eswitch = mean(curr_data(celltype.excit & celltype.pool_switch,:),1);
        Istay = mean(curr_data(celltype.inhib & celltype.pool_stay,:),1);
        Iswitch = mean(curr_data(celltype.inhib & celltype.pool_switch,:),1);
        
        %now find the crossover times
        all_groups = [Estay;Eswitch;Istay;Iswitch];%use legend labels as ref if you get confused
        Cpoints = find_crossover_points(all_groups);
        Yscale = cell2mat(Cpoints(:,2));
        Yscale = range(Yscale)*.05;
        hold on
        for Cidx = 1:numel(Cpoints(:,1))
            %my_color = tr_col{Cidx};
            Cp = Cpoints(Cidx,:); 
            %find the traces intersecting & which way they're going
            up_trace = cellfun(@(x) strcmp(x,'up'),Cp);
            up_trace = Cp{find(up_trace)-1};
            down_trace = cellfun(@(x) strcmp(x,'down'),Cp);
            down_trace = Cp{find(down_trace)-1};
            %plot in their respective colors 
            scatter(Cp{1},Cp{2}+Yscale,100,tr_col{up_trace},'^','Filled')
            scatter(Cp{1},Cp{2}-Yscale,100,tr_col{down_trace},'v','Filled')
        end
        
        l = plot(NaN,'Color',lcol,'LineWidth',lnsz); %for legend
        
        hold off
        Xlim_max = numel(curr_data(1,:));
        set(gca,'XLim',[0 Xlim_max]);
        Xticks = num2cell(get(gca,'Xtick'));
        %Xticks = Xticks(2:2:end-1);
        if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
        Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        legend(l,sprintf('%.0fHz',curr_net_info{j,3}),'location','best')
        
        if plt_idx == 9 || plt_idx == 10
            xlabel('Leave decision (ms)')
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',Nspeed),'Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('net #%i spiking (hz)',idx))
        end
          
    end
end
orient tall
int_fig_fn = [fig_fn '_xover'];

switch zoomed_fig
    case 'yes'
        %zoom in better
        Yl = arrayfun(@(x) x.Children,h,'UniformOutput',false);
        %skip legend & I-stay
        Yl = cellfun(@(x) x([2,4,5]),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) cat(1,x(:).YData),Yl,'UniformOutput',false);
        Yl = cellfun(@(x) max(x(:)),Yl);
        arrayfun(@(x,y) ylim(x,[0,y]),h,Yl,'UniformOutput',false);
        int_fig_fn = [int_fig_fn '_zoomed'];
    case 'no'
        axis tight
end

linkaxes(h,'x')
switch print_anything
    case 'yes'
        print(fullfile(fig_dir,int_fig_fn),'-djpeg','-r300')
end

close all;hold off



function cell_raster = sim_spikerate(cell_raster,timestep)

num_binsamps = 2e-3/timestep; %num samples in 2ms
raster_sz = size(cell_raster);
bin_magic = [numel(cell_raster(:,1)), num_binsamps, numel(cell_raster(1,:))/num_binsamps]; %set up for a magic trick
cell_raster = reshape(cell_raster,bin_magic);
cell_raster = squeeze(sum(cell_raster,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_raster = repmat(cell_raster,[num_binsamps 1 1]);
cell_raster = reshape(cell_raster,raster_sz); %put the rabit back in the hat
end




