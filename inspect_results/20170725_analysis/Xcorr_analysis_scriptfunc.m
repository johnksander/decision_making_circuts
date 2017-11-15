function  Xcorr_analysis_scriptfunc(celltype,targIDs,Xcor_method)
%NOTE: this is just a function version of Xcorr_analysis_20170729
%for making a bunch of figures...
%-----input:
% targIDs = {{'post-switch','stimulus'};{'pre-switch','stimulus'}};
% celltype = {'I-stay'};
% Xcor_method = 'coeff';

home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
addpath(fullfile(home_dir,'helper_functions'))
results_dir = fullfile(home_dir,'Results','change_switchspeed');
data_dir = fullfile(results_dir,'event_data');
output_dir = fullfile(results_dir,'Xcorr_analysis');
if strcmp(Xcor_method,'unbiased')
    output_dir = fullfile(output_dir,'unbiased_correction');
end
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

% %events2plot = {'stable','stable + stimulus','pre-switch (stimulus)','post-switch (stimulus)'};
% %now corresponding data IDs
% targIDs = {{'post-switch','stimulus'};{'pre-switch','stimulus'}};
% celltype = {'I-stay'};
% Xcor_method = 'coeff';

fontsz = 12;
mdl_colors = [[250 70 22]./255;0,0.4470,0.7410];
%Yax_labs = {'spike rate (Hz)';'per 2ms bin'};
%Xax_labs = 'state switch timecourse (ms)';
Yax_labs = '%s * %s';
Xax_labs = '%s lag offset (ms)';


%select data
cell_data = cell(numel(targIDs),num_models);
for select_idx = 1:numel(targIDs)
    ID2get = [targIDs{select_idx},celltype]; %just so func input below is a little easier...
    cell_data(select_idx,:) = ...
        grab_data(ID2get,dataIDs,model_names,event_labels,event_data);
end

%zero mean the spikerates
cell_data = cellfun(@(x) x - mean(x),cell_data,'UniformOutput',false);


orient landscape
for plot_idx = 1:num_models
    
    subplot(num_models,1,plot_idx)
    ln_color = mdl_colors(plot_idx,:);
    
    model_data = cell_data(:,plot_idx);
    
    %get cross correlation
    [r,laginds] = xcorr(model_data{1},model_data{2},Xcor_method); %do Xcorr
    
    plot(laginds,r,'Linewidth',2,'color',ln_color)
    
    set(gca,'XLim',[min(laginds)-1 max(laginds)+1]);
    Xticks = num2cell(get(gca,'Xtick'));
    Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
    set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    
    ylabel({sprintf(Yax_labs,targIDs{1}{1},targIDs{2}{1});...
        '(both with stimulus)'})
    xlabel(sprintf(Xax_labs,targIDs{2}{1}))
    
    title(model_names{plot_idx})
    set(gca,'Fontsize',fontsz)
end

%add a grand title at the top
grand_title([char(celltype) ' cells'],16)

%print it out
fig_XcorID = [strrep(targIDs{1}{1},'-','') 'X' strrep(targIDs{2}{1},'-','')];
fig_fn = strrep(char(celltype),'-','');
fig_fn = [fig_fn '_' fig_XcorID];

print(fullfile(output_dir,fig_fn),'-djpeg')
close all



    function output = grab_data(IDs,dataIDs,mdls,event_tags,data)
        
        output = cell(size(mdls));
        
        for mdl_idx = 1:numel(mdls)
            curr_mdl = mdls{mdl_idx};
            curr_event = cell2mat(cellfun(@(x) strcmp(x,IDs{1}),event_tags','UniformOutput',false));
            selection = dataIDs{curr_event}; %now get the right data for this event
            selection = cellfun(@(x) ... %cellfun setup
                [strcmp(x(:,1),curr_mdl),strcmp(x(:,2),IDs{2}),strcmp(x(:,3),IDs{3})],...%function
                selection,'UniformOutput',false);
            selection = cellfun(@(x,y) x == y,...
                cellfun(@(x) sum(x,2),selection,'UniformOutput',false),...
                cellfun(@(x) repmat(size(x,2),size(x,1),1),selection,'UniformOutput',false),...
                'UniformOutput',false); %major nested cellfun alert
            selection = cellfun(@(x,y) x(y,:),data{curr_event},selection,'UniformOutput',false);
            selection = cell2mat(selection); %got you!!
            output{mdl_idx} = selection;
        end
        
    end


    function grand_title(txt,sz)
        T=axes('Units','Normal','Position',[-.2 .09 .85 .85],'Visible','off');
        set(get(T,'Title'),'Visible','on')
        title(txt,'FontSize',sz,'FontWeight','b');
    end


end