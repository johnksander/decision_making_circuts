clear
clc
format compact
hold off;close all

rescale_plane = 'on';
outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu'

%result summaries
fontsz = 16;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_baseline';
basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
%figdir = fullfile(basedir,'Results',[sim_name '_figures']);
figdir = fullfile(basedir,'Results','parsweep_stims_figures');
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
num_files = numel(output_fns);

%network parameter info 
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
network_pair_info = load(fullfile(basedir,'helper_functions','network_pairs'));
NPI_durations = network_pair_info.durations;
network_pair_info = network_pair_info.network_pairs;
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],...
    network_pair_info,'UniformOutput',false);
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],...
    network_pair_info,'UniformOutput',false);
network_pair_info = cellfun(@(x,y) [x,y],network_pair_info,NPI_durations,'UniformOutput',false);

%get results
file_data = cell(num_files,2);
for idx = 1:num_files
    curr_file = load(output_fns{idx});
    %store parameters
    file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{:};
    %just get all of them, baseline test
    state_durations = cellfun(@(x) x(:,1),state_durations,'UniformOutput',false);
    state_durations = vertcat(state_durations{:});
    state_durations = vertcat(state_durations{:}); %ooo that's annoying
    %convert to time
    state_durations = state_durations * timestep;
    file_data{idx,1} = state_durations;
end

%search for jobs with identical parameters, collapse distributions
%get the randomized network parameters
job_params = cellfun(@(x)  [x.ItoE, x.EtoI],file_data(:,2),'UniformOutput',false);
job_params = vertcat(job_params{:});
uniq_params = unique(job_params,'rows');
num_jobs = size(uniq_params,1);
fprintf('----------------------\n')
fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)

%collapse duplicate job parameters
result_data = cell(num_jobs,2);
for idx = 1:num_jobs
    %find all matching
    curr_file = ismember(job_params,uniq_params(idx,:),'rows');
    %collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
end



%just grab some simple stats here
num_states = cellfun(@(x) numel(x),result_data(:,1));
mu_duration = cellfun(@(x) mean(x),result_data(:,1));
med_duration = cellfun(@(x) median(x),result_data(:,1));
logmu_dur = cellfun(@(x) mean(log(x)),result_data(:,1));
%get the network parameters
EtoE = cellfun(@(x)  x.EtoE,result_data(:,2));
ItoE = cellfun(@(x)  x.ItoE,result_data(:,2));
EtoI = cellfun(@(x)  x.EtoI,result_data(:,2));


switch outcome_stat
    case 'mu'
        outcome = mu_duration;
        Zlabel = 'mean duration (s)';
        figdir = fullfile(figdir,'mean_duration');
    case 'med'
        outcome = med_duration;
        Zlabel = 'median duration (s)';
        figdir = fullfile(figdir,'med_duration');
    case 'logmu'
        network_pair_info = cellfun(@(x) [x(:,1:4) num2cell(log(cell2mat(x(:,end))))],...
            network_pair_info,'UniformOutput',false);
        outcome = logmu_dur;
        Zlabel = {'mean stay duration';'log(s)'};
        figdir = fullfile(figdir,'logmean_duration');
end

if ~isdir(figdir),mkdir(figdir);end



hold off;close all
%try countour plot
xlin = linspace(min(ItoE),max(ItoE),num_jobs);
ylin = linspace(min(EtoI),max(EtoI),num_jobs);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(ItoE,EtoI,outcome);
f.Method = 'natural';
Z = f(X,Y);
switch rescale_plane
    case 'on'
        Z(Z < 0) = 0;
        Z(Z > max(outcome)) = max(outcome);
end
mask = tril(ones(size(Z)),100);
mask = logical(flipud(mask));
Z(~mask) = NaN;
figure
[c,h] = contour(X,Y,Z,'LineWidth',2);
levs = h.LevelList;
h.LevelList = levs(1:2:end);
clabel(c,h);
title('State durations')
xlabel('I-to-E strength')
ylabel('E-to-I strength')
set(gca,'FontSize',fontsz)
hold on
pointsz = 200;
for idx = 1:numel(network_pair_info)
    curr_network = network_pair_info{idx};
    for Nidx = 1:2
       scatter(curr_network{Nidx,1},curr_network{Nidx,2},pointsz,'red','filled','p')
    end
    adj_y = (max(ylin) -  curr_network{2,2}) * .125;
    ajd_x = (max(xlin) -  curr_network{2,1}) * .15;
    text(curr_network{2,1}+ajd_x,curr_network{2,2}+adj_y,sprintf('Net #%i',idx),'color','red',...
        'HorizontalAlignment','center','Fontweight','b','Fontsize',fontsz)
    
end
%print(fullfile(figdir,'contour_marked'),'-djpeg')


Nplots = cell2mat(cellfun(@(x) cell2mat(x(:,[1,2,5])),network_pair_info,'UniformOutput',false));



%this 3D contour plot is really good actually
hold off;close all
H = figure;
[c,h] = contour3(X,Y,Z,'LineWidth',5);
levs = h.LevelList;
h.LevelList = levs(1:2:end);
%clabel(c,h);
axis tight; hold on
plot3(ItoE,EtoI,outcome,'.','MarkerSize',15,'color','blue') %give them outlines
plot3(ItoE,EtoI,outcome,'.','MarkerSize',10,'color','green') 
plot3(Nplots(:,1),Nplots(:,2),Nplots(:,3),'p','MarkerSize',25,'MarkerFaceColor','red') 
xlabel('I-to-E strength')
ylabel('E-to-I strength')
zlabel(Zlabel)
set(gca,'FontSize',fontsz)
%savefig(H,fullfile(figdir,'contour3D_plot_marked'),'compact')






%surface plot
hold off;close all
fontsz = 20;
plt_alph = .75;
num_gridlines = 400;
xlin = linspace(min(ItoE),max(ItoE),num_gridlines);
ylin = linspace(min(EtoI),max(EtoI),num_gridlines);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(ItoE,EtoI,outcome);
f.Method = 'natural';
Z = f(X,Y);
mask = tril(ones(size(Z)),30);
mask = logical(flipud(mask));
Z(~mask) = NaN;
figure;
mesh(X,Y,Z,'FaceAlpha',plt_alph,'EdgeAlpha',plt_alph) %interpolated
axis tight; hold on

scatter3(ItoE,EtoI,outcome,30,'blue','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)  %give them outlines
scatter3(ItoE,EtoI,outcome,15,'green','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1) 
scatter3(Nplots(:,1),Nplots(:,2),Nplots(:,3),1000,'red','p','filled',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1); %mark 'em 
%plot3(ItoE,EtoI,outcome,'.','MarkerSize',15,'color','blue') %give them outlines
%plot3(ItoE,EtoI,outcome,'.','MarkerSize',10,'color','green') 
%plot3(Nplots(:,1),Nplots(:,2),Nplots(:,3),'p','MarkerSize',25,'MarkerFaceColor','red') 
xlabel({'within-pool inhibition';'(I-to-E strength)'},'FontWeight','b')
ylabel({'cross-pool inhibition';'(E-to-I strength)'},'FontWeight','b')
zlabel(Zlabel,'FontWeight','b')
set(gca,'FontSize',fontsz)
view(-45,27)
hidden off
rotate3d on
%also see--
%https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html#bsovi2t
switch rescale_plane
    case 'on'
        Zlim = get(gca,'ZLim');
        Zlim(1) = 0;
        Zlim(2) = max(outcome);
        set(gca,'ZLim',Zlim);
        caxis(Zlim)
end
%for brownbag
%make it big
set(gcf,'units','normalized','outerposition',[0 0 .75 1])
%fix the labels 
xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
xh_pos = get(xh, 'Position');
set(xh, 'Position',xh_pos+[.05,.1,0],'Rotation',20.5)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
yh_pos = get(yh, 'Position');
set(yh, 'Position',yh_pos+ [-.075,.15,0],'Rotation',-19)
%print high-res
print(fullfile(figdir,'brownbag_marked_surface_plot'),'-djpeg','-r300')
%https://www.mathworks.com/matlabcentral/answers/41800-remove-sidewalls-from-surface-plots
%set(gca,'FontSize',fontsz)
%savefig(fullfile(figdir,'surface_plot_marked'))

keyboard