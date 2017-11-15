clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.
%1) this code takes datafiles made from inspect_results_event_data

home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
addpath(fullfile(home_dir,'helper_functions'))
results_dir = fullfile(home_dir,'Results','change_switchspeed');
data_dir = fullfile(results_dir,'event_data');
output_dir = fullfile(results_dir,'network_behavior_figs');
if ~isdir(output_dir), mkdir(output_dir); end

%load in preprocessed data
dataset = load(fullfile(data_dir,'event_data.mat'));
event_data = dataset.event_data;
event_windows = dataset.event_windows;
adjusted_switch_onsets = dataset.adjusted_switch_onsets;
cell_inds = dataset.info.cell_inds;
timestep = dataset.info.timestep;
recorded_switchtime = dataset.info.recorded_switchtime;
event_labels = dataset.info.event_labels;
sim_names = dataset.info.sim_names;




%LOOK AT cell_inds FOR CORRECT ORDERING OF THESE LABELS
cell_inds;
cell_labels = {'E-stay','I-stay','E-switch','I-switch'};
cell_inds = cell2mat(struct2cell(cell_inds));
clear dataset


%do plots
model_names = cellfun(@(x) strsplit(x,'_'), sim_names,'UniformOutput',false);
model_names = cellfun(@(x) x{1}, model_names,'UniformOutput',false);
model_names = unique(model_names);

num_events = numel(event_data);
num_models = numel(model_names);
num_celltypes = numel(cell_labels);

dataIDs = make_dataIDs(event_data,sim_names,cell_labels); %make ID structure matching event_data

events2plot = {'stable','stable + stimulus','pre-switch (stimulus)','post-switch (stimulus)'};
%now corresponding data IDs
plotIDs = {{'stable','baseline'};{'stable','stimulus'};{'pre-switch','stimulus'};{'post-switch','stimulus'}};

fontsz = 12;
mdl_colors = [[250 70 22]./255;0,0.4470,0.7410];
Yax_labs = {'spike rate (Hz)';'per 2ms bin'};
Xax_labs = 'state switch timecourse (ms)';

for figidx = 1:num_celltypes
    
    cell2plot = cell_labels{figidx};
    num_subplots = num_models * numel(events2plot);
    
    orient landscape
    for plot_idx = 1:num_subplots
        
        plot_row_idx = floor((plot_idx + 1)/num_models);
        
        if mod(plot_idx,2) == 1
            mdl_idx = 1; %I just wanted it in this order... annoying without an if-statement
        elseif mod(plot_idx,2) == 0
            mdl_idx = 2;
        end
        
        ax(plot_idx) = subplot(numel(events2plot),num_models,plot_idx);
        
        model2plot = model_names{mdl_idx}; %I just wanted it in this order... annoying without an if-statement
        ln_color = mdl_colors(mdl_idx,:);
        plot_type = plotIDs{plot_row_idx};
        curr_event = cell2mat(cellfun(@(x) strcmp(x,plot_type{1}),event_labels','UniformOutput',false));
        %find the correct data for the current plot
        timecourse_data = dataIDs{curr_event};
        timecourse_data = cellfun(@(x) ... %cellfun setup
            [strcmp(x(:,1),model2plot),strcmp(x(:,2),plot_type{2}),strcmp(x(:,3),cell2plot)],...%function
            timecourse_data,'UniformOutput',false);
        timecourse_data = cellfun(@(x,y) x == y,...
            cellfun(@(x) sum(x,2),timecourse_data,'UniformOutput',false),...
            cellfun(@(x) repmat(size(x,2),size(x,1),1),timecourse_data,'UniformOutput',false),...
            'UniformOutput',false); %major nested cellfun alert
        timecourse_data = cellfun(@(x,y) x(y,:),event_data{curr_event},timecourse_data,'UniformOutput',false);
        timecourse_data = cell2mat(timecourse_data); %got you!!
        
        plot(timecourse_data,'Linewidth',2,'color',ln_color)
        
        %Xlim_max = numel(timecourse_data(1,:));
        %set(gca,'XLim',[0 Xlim_max]);
        Xticks = num2cell(get(gca,'Xtick'));
        if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
        adj_switch = adjusted_switch_onsets(curr_event);
        Xlabs = cellfun(@(x) sprintf('%+i',((x-adj_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        if mod(plot_row_idx,numel(events2plot)) == 0
            xlabel(Xax_labs)
        end
        if mod(plot_idx,num_models) == 1
            ylabel(Yax_labs)
        end
        if plot_row_idx == 1
            place_legend(cell2plot,model2plot)
        end
        
        title(events2plot{plot_row_idx})
        set(gca,'Fontsize',fontsz)
        
    end
    
    %add a grand title at the top
    grand_title([cell2plot ' cells'],16)
    
    %align all the axes within models
    linkaxes(ax(1:2:end),'y');
    linkaxes(ax(2:2:end),'y');
    
    %print it out
    fig_fn = strrep(cell2plot,'-','_');
    print(fullfile(output_dir,fig_fn),'-djpeg')
    close all
end

function place_legend(cell2plot,model2plot)

%this is annoying...
if strcmp(cell2plot,'E-switch') | strcmp(cell2plot,'I-stay')
    legend(model2plot,'location','northeast') %these are fine
elseif strcmp(cell2plot,'E-stay') | strcmp(cell2plot,'I-switch')
    axpos = get(gca,'Position');
    leg = legend(model2plot,'location','southeast');
    legpos = get(leg,'Position'); %keep width & height
    legpos(1) = (axpos(1) + axpos(3)) - legpos(3) - .02;
    legpos(2) = axpos(2) + .02;
    set(leg,'Position',legpos,'Units','normalized');
end
end


function grand_title(txt,sz)
T=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(T,'Title'),'Visible','on')
title(txt,'FontSize',sz,'FontWeight','b');
end