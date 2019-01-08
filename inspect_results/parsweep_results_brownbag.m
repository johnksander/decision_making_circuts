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
figdir = fullfile(basedir,'Results',[sim_name '_figures']);
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
num_files = numel(output_fns);
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
        outcome = logmu_dur;
        Zlabel = {'mean stay duration';'log(s)'};
        figdir = fullfile(figdir,'logmean_duration');
end

if ~isdir(figdir),mkdir(figdir);end


%inhibition scatter
net_inhibition = ItoE.*EtoI;
%EtoI & ItoE parameter ranges are different (EtoI goes lower, ItoE goes much higher)
%plotting so that the ranges are equal
NI_range = EtoI >= min(ItoE) & ItoE <= max(EtoI);
h1 = subplot(2,2,1);
larger_param = ItoE > EtoI & NI_range;
scatter(net_inhibition(larger_param),outcome(larger_param));
title('I-to-E > E-to-I')
ylabel(Zlabel)
xlabel('net inhibition')
set(gca,'FontSize',12)
h2 = subplot(2,2,3);
larger_param = ItoE < EtoI & NI_range;
scatter(net_inhibition(larger_param),outcome(larger_param));
title('I-to-E < E-to-I')
ylabel(Zlabel)
xlabel('net inhibition')
set(gca,'FontSize',12)
linkaxes([h1,h2],'xy')
axis tight

h1 = subplot(2,2,2);
larger_param = ItoE > EtoI & NI_range;
histogram(ItoE(larger_param));
hold on
histogram(EtoI(larger_param));
hold off
legend('I-to-E','E-to-I')
title('parameter range')
ylabel('frequency')
xlabel('strength')
set(gca,'FontSize',12)
h2 = subplot(2,2,4);
larger_param = ItoE < EtoI & NI_range;
histogram(ItoE(larger_param));
hold on
histogram(EtoI(larger_param));
hold off
legend('I-to-E','E-to-I')
title('parameter range')
ylabel('frequency')
xlabel('strength')
set(gca,'FontSize',12)
linkaxes([h1,h2],'xy')
axis tight
%print(fullfile(figdir,'net_inhib_scatter'),'-djpeg')



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
mesh(X,Y,Z,'FaceAlpha',plt_alph,'EdgeAlpha',plt_alph) %interpolated
axis tight; hold on
%these feel like they need to be slightly bigger 
%scatter3(ItoE,EtoI,outcome,25,'black','filled','MarkerEdgeAlpha',1,'MarkerEdgeAlpha',1) %give them outlines
scatter3(ItoE,EtoI,outcome,30,'red','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1) 
xlabel({'within-pool inhibition';'(I-to-E strength)'},'FontWeight','b')
ylabel({'cross-pool inhibition';'(E-to-I strength)'},'FontWeight','b')
zlabel(Zlabel,'FontWeight','b')
set(gca,'FontSize',fontsz)
view(-45,27)
%view(-24,24)
%setting view to overhead gives heatmap!!
%view(0,90) %that looks good too
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
%savefig(fullfile(figdir,'surface_plot'))
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
print(fullfile(figdir,'brownbag_surface_plot'),'-djpeg','-r300')

%https://www.mathworks.com/matlabcentral/answers/41800-remove-sidewalls-from-surface-plots
%savefig(fullfile(figdir,'surface_plot'))
%view(-24,24)

figure
hold off

xlin = linspace(min(ItoE),max(ItoE),num_jobs);
ylin = linspace(min(EtoI),max(EtoI),num_jobs);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(ItoE,EtoI,outcome);
f.Method = 'natural';
Z = f(X,Y);
Z = flipud(Z); %so the bottom left corner is low x, low y
switch rescale_plane
    case 'on'
        Z(Z < 0) = 0;
        Z(Z > max(outcome)) = max(outcome);
end
mask = tril(ones(size(Z)),100);
mask = logical(mask);
Z(~mask) = NaN;
colormap(parula)
imagesc(Z)
xlabel('I-to-E strength')
ylabel('E-to-I strength')
title('State durations')
set(gca,'Fontsize',fontsz)
Xticks = get(gca,'Xtick');
Xticks = linspace(min(ItoE),max(ItoE),numel(Xticks));
Xlabs = cellfun(@(x) sprintf('%.2f',x),num2cell(Xticks),'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs)
Yticks = get(gca,'Ytick');
Yticks = linspace(max(EtoI),min(EtoI),numel(Yticks));
Ylabs = cellfun(@(x) sprintf('%.2f',x),num2cell(Yticks),'UniformOutput', false);
set(gca, 'YTickLabel', Ylabs')
colb = colorbar;
colb.Label.String = Zlabel;
colb.FontSize = 16;
set(gca,'FontSize',fontsz)
print(fullfile(figdir,'heatmap'),'-djpeg')


% %this is okay.....
% sz = linspace(25,400,num_jobs);
% c = colormap(summer(num_jobs));
% [~,inds] = sort(outcome,'descend');
% c = c(inds,:);
% sz = sz(inds);
% figure
% scatter(ItoE,EtoI,sz,c,'filled')
% colb = colorbar;
% Tlabs = colb.TickLabels;
% Tlabs = linspace(min(outcome),max(outcome),numel(Tlabs));
% Tlabs = num2cell(Tlabs);
% Tlabs = cellfun(@(x) sprintf('%.2f',x),Tlabs,'UniformOutput', false);
% colb.TickLabels = Tlabs;
% xlabel('I-to-E strength')
% ylabel('E-to-I strength')


hold off


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
print(fullfile(figdir,'contour'),'-djpeg')


%this 3D contour plot is really good actually
H = figure;
[c,h] = contour3(X,Y,Z,'LineWidth',5);
levs = h.LevelList;
h.LevelList = levs(1:2:end);
%clabel(c,h);
axis tight; hold on
plot3(ItoE,EtoI,outcome,'.','MarkerSize',25,'color','red') %nonuniform
xlabel('I-to-E strength')
ylabel('E-to-I strength')
zlabel(Zlabel)
set(gca,'FontSize',fontsz)
savefig(H,fullfile(figdir,'contour3D_plot'),'compact')


hold off
figure
histogram(outcome)
ylabel('frequency')
xlabel(Zlabel)
title('State durations')
set(gca,'FontSize',fontsz)
print(fullfile(figdir,'histogram'),'-djpeg')

% 
% searchlight_radius = .05;
% mesh2keep = zeros(size(X));
% for idx = 1:num_jobs
%     x = ItoE(idx);
%     y = EtoI(idx);
%     z = outcome(idx);
% %     sphere_voxels = logical((Y - y(1)).^2 + ...
% %         (X - x(1)).^2 + (Z - z(1)).^2 ...
% %         <= searchlight_radius.^2); %adds a logical
% 
% sphere_voxels = (Y - y(1)) <= .01 & ...
%     (X - x(1)) <= .01 & (Z - z(1)) <= .1;
% 
%     mesh2keep(sphere_voxels) = 1;
% end
% %imagesc(mesh2keep)
% 
% X(~mesh2keep) = NaN;
% Y(~mesh2keep) = NaN;
% Z(~mesh2keep) = NaN;
% 

%this does grid just around the dataponts, doesn't look great tho..

% switch do_surplot
%     case 'on'
%         hold off;close all
%
%         figure(1)
%         xlin = linspace(min(ItoE),max(ItoE),num_files);
%         ylin = linspace(min(EtoI),max(EtoI),num_files);
%         [X,Y] = meshgrid(xlin,ylin);
%         Z = griddata(ItoE,EtoI,mu_duration,X,Y);
%         figure
%         mesh(X,Y,Z) %interpolated
%         axis tight; hold on
%         plot3(ItoE,EtoI,mu_duration,'.','MarkerSize',25,'color','red') %nonuniform
%         xlabel('ItoE')
%         ylabel('EtoI')
%         zlabel('mu duration')
%         %setting view to overhead gives heatmap!!
%         %view(0,90) %that looks good too
%         hidden off
%         rotate3d on
%
%         %also see--
%         %https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html#bsovi2t
%         switch rescale_plane
%             case 'on'
%                 Zlim = get(gca,'ZLim');
%                 Zlim(1) = 0;
%                 Zlim(2) = max(mu_duration);
%                 set(gca,'ZLim',Zlim);
%                 caxis(Zlim)
%         end
% end
