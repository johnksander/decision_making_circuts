clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.
%1) this code takes datafiles made from inspect_results_event_data

home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
addpath(fullfile(home_dir,'RDMfuncs'))
results_dir = fullfile(home_dir,'Results','change_switchspeed');
data_dir = fullfile(results_dir,'event_data');
output_dir = fullfile(results_dir,'RSA');
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
%cell_labels = {'E-stay','I-stay','E-switch','I-switch'};
cell_labels = {'E-st','I-st','E-sw','I-sw'}; %for plotting
clear dataset

model_names = cellfun(@(x) strsplit(x,'_'), sim_names,'UniformOutput',false);
model_names = cellfun(@(x) x{1}, model_names,'UniformOutput',false);
model_names = unique(model_names);
num_models = numel(model_names);

num_events = numel(event_data);
num_celltypes = numel(fieldnames(cell_inds));

num_entries = num_events * num_celltypes * 2; %a little hardcoded, there's two model runs per event

RDMs = NaN(num_entries,num_entries,num_models);

for idx = 1:num_models
    
    model_inds = strncmpi(model_names{idx},sim_names,numel(model_names{idx}));
    model_data = cellfun(@(x) x(model_inds),event_data,'UniformOutput',false);
    model_data = cat(1,model_data{:});
    model_data = cat(1,model_data{:});
    
    RDMs(:,:,idx) = corr(model_data','type','kendall');
end


vec_mask = logical(triu(ones(size(RDMs(:,:,1))),1));

RDMvecs = NaN(sum(vec_mask(:)),num_models);
for idx = 1:num_models
    currRDM = RDMs(:,:,idx);
    RDMvecs(:,idx) = currRDM(vec_mask);
end

true_fit =  corr(RDMvecs,'type','kendall','rows','pairwise');
true_fit = unique(true_fit(~logical(eye(size(true_fit))))); %lol
disp(sprintf('Kendall''s tau = %.4f',true_fit))


%do permutations, harcoding for 1 vs 1 RDM right now...
%WAIT... what're you doing man. You don't need to test if they're related...
% num_perms = 10000;
% RDM2permute = RDMs(:,:,1);
% testRDM = RDMvecs(:,2);
%
% permutedRDMs = NaN(sum(vec_mask(:)),num_perms);
% for idx = 1:num_perms
%     curr_order = randperm(num_entries);
%     currRDM = RDM2permute(curr_order,curr_order);
%     permutedRDMs(:,idx) = currRDM(vec_mask);
% end
%
% null_dist = corr(permutedRDMs,testRDM,'type','kendall','rows','pairwise');
% true_fit =  corr(RDMvecs,'type','kendall','rows','pairwise');
% true_fit = unique(true_fit(~logical(eye(size(true_fit))))); %lol



%visualization time
%super hardcoded now need to check cell_labels var for correctness!!!!
cell_labels;

txtsz = 10;
close all

for idx = 1:num_models
    
    currRDM =  RDMs(:,:,idx); %make visualization code somewhat easier...
    alpha = ~isnan(currRDM);
    fig = image(scale01(rankTransform_equalsStayEqual(currRDM,1)),'CDataMapping','scaled','AlphaData',alpha);
    set(gca,'CLim',[0 1],'CLimMode','manual');
    colormap('parula')
    cbar = colorbar;
    cbar.Label.String = 'entry rank';
    cbar.FontSize = 16;
    
    figpos = get(gca,'Position');
    figpos(2) = figpos(2) + .0225;
    set(gca,'Position',figpos);
    
    XLims = get(gca,'XLim');
    YLims = get(gca,'YLim');
    label_inds = linspace(figpos(2),figpos(4),9);
    label_inds = label_inds(2:3:end);
    
    title([model_names{idx} ' model'],'FontSize',16,'FontWeight','bold')
    
    XTick = 1:num_entries;
    set(gca,'XTick',XTick);
    YTick = 1:num_entries;
    set(gca,'YTick',YTick);
    
    cell_type_label = num2cell(YTick)'; %again, check cell_labels var for correctness
    odds = logical(cellfun(@(x) mod(x,2),cell_type_label));
    cell_type_label(~odds) = {'In'};
    cell_type_label(odds) = {'Ex'};
    
    inds4stay = 1:2:num_entries;
    inds4stay = inds4stay(1:2:end);
    inds4switch = 1:2:num_entries;
    inds4switch = inds4switch(2:2:end);
    cell_pool_label = cell(size(cell_type_label));
    cell_pool_label(inds4stay) = {'stay'};
    cell_pool_label(inds4switch) = {'switch'};
    cell_pool_label(2:2:end) = {''};
    cell_pool_label = char(pad(cell_pool_label,'left'));
    
    inds4event = 4:4:num_entries;
    inds4event = inds4event(1:2:end);
    fig_event_label = cell(size(cell_type_label));
    fig_event_label(inds4event) = event_labels;
    fig_event_label(cellfun(@isempty,fig_event_label)) = {''};
    fig_event_label = char(pad(fig_event_label,'Both'));
    
    inds4stim = zeros(num_entries,1);
    inds4stim(5:8:end) = 1;%this is so dumb...
    inds4stim(6:8:end) = 1;%this is so dumb...
    inds4stim(7:8:end) = 1;%this is so dumb...
    inds4stim(8:8:end) = 1;%this is so dumb...
    inds4stim = logical(inds4stim);
    
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    offset = 0.0225;
    x = XLims(1) - offset*diff(XLims);
    y = YLims(2) + offset*diff(YLims);
    Ytext = text(repmat(x,size(YTick)),YTick,cat(1,cell_type_label{:}),'HorizontalAlignment','Center',...
        'FontSize',txtsz,'Rotation',90);
    Xtext = text(XTick,repmat(y,size(XTick)),cat(1,cell_type_label{:}),'HorizontalAlignment','Center',...
        'FontSize',txtsz);
    
    paint_it_red(Ytext,inds4stim);
    paint_it_red(Xtext,inds4stim);
    
    %offset = 0.07;
    x = XLims(1) - .05*diff(XLims);
    y = YLims(2) + .11*diff(YLims);
    Ytext = text(repmat(x,size(YTick)),YTick,cell_pool_label,'HorizontalAlignment','Right',...
        'FontSize',txtsz,'Rotation',45);
    Xtext = text(XTick,repmat(y,size(XTick)),cell_pool_label,'HorizontalAlignment','Left',...
        'FontSize',txtsz,'Rotation',45);
    
    paint_it_red(Ytext,inds4stim);
    paint_it_red(Xtext,inds4stim);
    
    offset = .135;
    x = XLims(1) - offset*diff(XLims);
    y = YLims(2) + offset*diff(YLims);
    Ytext = text(repmat(x,size(YTick)),YTick +.5,fig_event_label,'HorizontalAlignment','Center',...
        'FontSize',txtsz,'Rotation',90);
    Xtext = text(XTick + .5,repmat(y,size(XTick)),fig_event_label,'HorizontalAlignment','Center',...
        'FontSize',txtsz);
    
    paint_it_red(Ytext,inds4stim);
    paint_it_red(Xtext,inds4stim);
    
    print(fullfile(output_dir,[model_names{idx} '_model']),'-djpeg')
    close all
end

function paint_it_red(h,i)
h = h(i);
set(h,'Color','r');
end


%
% cell_row = cellfun(@(x) [' ' x],cell_labels,'UniformOutput',false);
% cell_row{1} = strrep(cell_row{1},' ',''); %even out the other side
% cell_row = cellfun(@(x) strrep(x,' ',' | '),cell_row,'UniformOutput',false);
% cell_row = cat(2,cell_row{:});
%
% run_type_row = 'baseline | stimulus'; %dumb and hardcoded!
%
% label_inds = linspace(0,1,9);
% label_inds = label_inds(2:3:end);
%
% txtsz = 10;
%
% for idx = 1:num_models
%
%     currRDM =  RDMs(:,:,idx); %make visualization code somewhat easier...
%     alpha = ~isnan(currRDM);
%     image(scale01(rankTransform_equalsStayEqual(currRDM,1)),'CDataMapping','scaled','AlphaData',alpha);
%     set(gca,'CLim',[0 1],'CLimMode','manual');
%     colormap('parula')
%     cbar = colorbar;
%     cbar.Label.String = 'entry rank';
%     set(gca,'XTick',[]);
%     set(gca,'YTick',[]);
%     for label_idx = 1:numel(label_inds)
%         text(label_inds(label_idx),-0.03,cell_row,'Units','normalized','HorizontalAlignment','center',...
%             'FontSize',txtsz)
%         text(label_inds(label_idx),-0.07,run_type_row,'Units','normalized','HorizontalAlignment','center',...
%             'FontSize',txtsz)
%         text(label_inds(label_idx),-0.11,event_labels{label_idx},'Units','normalized',...
%             'HorizontalAlignment','center','FontSize',txtsz)
%         text(-0.03,label_inds(label_idx),cell_row,'Units','normalized','HorizontalAlignment','center',...
%             'Rotation',90,'FontSize',txtsz)
%         text(-0.07,label_inds(label_idx),run_type_row,'Units','normalized',...
%             'HorizontalAlignment','center','Rotation',90,'FontSize',txtsz)
%         text(-0.11,label_inds(label_idx),event_labels{label_idx},'Units','normalized',...
%             'HorizontalAlignment','center','Rotation',90,'FontSize',txtsz)
%     end
%     title([model_names{idx} ' model'],'FontWeight','bold')
%
%     set(gca,'FontSize',16);
%
%     %print(fullfile(output_dir,[model_names{idx} '_model']),'-djpeg')
%
% end
%



%
% %super hardcoded now need to check cell_labels var for correctness!!!!
% cell_type_label = 'Ex  |  In  |  Ex  |  In';
% cell_pool_label =    'stay    |    switch';
% run_type_row = 'baseline  |  stimulus'; %dumb and hardcoded!
%
% label_inds = linspace(0,1,9);
% label_inds = label_inds(2:3:end);
%
% txtsz = 10;
% close all
% for idx = 1:num_models
%
%     currRDM =  RDMs(:,:,idx); %make visualization code somewhat easier...
%     alpha = ~isnan(currRDM);
%     image(scale01(rankTransform_equalsStayEqual(currRDM,1)),'CDataMapping','scaled','AlphaData',alpha);
%     set(gca,'CLim',[0 1],'CLimMode','manual');
%     colormap('parula')
%     cbar = colorbar;
%     cbar.Label.String = 'entry rank';
%     set(gca,'XTick',[]);
%     set(gca,'YTick',[]);
%     figpos = get(gca,'Position');
%     figpos(2) = figpos(2) + .0225;
%     set(gca,'Position',figpos);
%
%     for label_idx = 1:numel(label_inds)
%         text(label_inds(label_idx),-0.02,cell_type_label,'Units','normalized','HorizontalAlignment','center',...
%             'FontSize',txtsz)
%         text(label_inds(label_idx),-0.06,cell_pool_label,'Units','normalized','HorizontalAlignment','center',...
%             'FontSize',txtsz)
%         text(label_inds(label_idx),-0.10,run_type_row,'Units','normalized','HorizontalAlignment','center',...
%             'FontSize',txtsz)
%         text(label_inds(label_idx),-0.14,event_labels{label_idx},'Units','normalized',...
%             'HorizontalAlignment','center','FontSize',txtsz)
%
%         text(-0.02,label_inds(label_idx),cell_type_label,'Units','normalized','HorizontalAlignment','center',...
%             'Rotation',90,'FontSize',txtsz)
%         text(-0.06,label_inds(label_idx),cell_pool_label,'Units','normalized','HorizontalAlignment','center',...
%             'Rotation',90,'FontSize',txtsz)
%         text(-0.10,label_inds(label_idx),run_type_row,'Units','normalized',...
%             'HorizontalAlignment','center','Rotation',90,'FontSize',txtsz)
%         text(-0.14,label_inds(label_idx),event_labels{label_idx},'Units','normalized',...
%             'HorizontalAlignment','center','Rotation',90,'FontSize',txtsz)
%     end
%     title([model_names{idx} ' model'],'FontWeight','bold')
%
%     set(gca,'FontSize',16);
%
%     %print(fullfile(output_dir,['newlabel_' model_names{idx} '_model']),'-djpeg')
%     close all
% end
