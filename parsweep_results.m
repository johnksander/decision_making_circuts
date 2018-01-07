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
job_params = cellfun(@(x)  [x.ItoE, x.ItoE],file_data(:,2),'UniformOutput',false);
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
    case 'med'
        outcome = med_duration;
    case 'logmu'
        outcome = logmu_dur;
end


%surface plot

xlin = linspace(min(ItoE),max(ItoE),400);
ylin = linspace(min(EtoI),max(EtoI),400);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(ItoE,EtoI,outcome);
f.Method = 'natural';
Z = f(X,Y);
mask = tril(ones(size(Z)),30);
mask = logical(flipud(mask));
Z(~mask) = NaN;
figure
mesh(X,Y,Z) %interpolated
axis tight; hold on
plot3(ItoE,EtoI,outcome,'.','MarkerSize',25,'color','red') %nonuniform
xlabel('ItoE')
ylabel('EtoI')
zlabel('mu duration')
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
%https://www.mathworks.com/matlabcentral/answers/41800-remove-sidewalls-from-surface-plots


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
figure
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
colb.Label.String = 'mean duration (s)';
colb.FontSize = 16;


hold off


%this is okay.....
sz = linspace(25,400,num_jobs);
c = colormap(summer(num_jobs));
[~,inds] = sort(outcome,'descend');
c = c(inds,:);
sz = sz(inds);
figure
scatter(ItoE,EtoI,sz,c,'filled')
colb = colorbar;
Tlabs = colb.TickLabels;
Tlabs = linspace(min(outcome),max(outcome),numel(Tlabs));
Tlabs = num2cell(Tlabs);
Tlabs = cellfun(@(x) sprintf('%.2f',x),Tlabs,'UniformOutput', false);
colb.TickLabels = Tlabs;
xlabel('I-to-E strength')
ylabel('E-to-I strength')


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

xlabel('I-to-E strength')
ylabel('E-to-I strength')

%this 3D contour plot is really good actually
figure
[c,h] = contour3(X,Y,Z,'LineWidth',5);
levs = h.LevelList;
h.LevelList = levs(1:2:end);
%clabel(c,h);
axis tight; hold on
plot3(ItoE,EtoI,outcome,'.','MarkerSize',25,'color','red') %nonuniform

xlabel('I-to-E strength')
ylabel('E-to-I strength')
zlabel('mu duration')



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
