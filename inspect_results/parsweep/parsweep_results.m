clear
clc
format compact
hold off;close all

num_workers = 24; %for parfor loading results 
rescale_plane = 'on';
outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu' ||| 'E-rate' | 'I-rate'

%result summaries
fontsz = 16;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..
limit_prange = 'no'; %'yes' | 'no' if yes, set maximums for connection parameters

EI_max = .75; IE_max = 8; %these set connection maximums
%for the fastD
%EI_max = .35; IE_max = 3.75; %these set connection maximums

%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_fastD_Rlim_baseline';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
figdir = fullfile(basedir,'Results',['figures_' sim_name]);
resdir = fullfile(basedir,'Results',sim_name);
addpath(fullfile(basedir,'helper_functions'))
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
%get results
num_files = numel(output_fns);
file_data = cell(num_files,3);
%parfor stuff
output_log = fullfile(resdir,'output_log.txt');
special_progress_tracker = fullfile(resdir,'SPT.txt');
if exist(special_progress_tracker) > 0
    delete(special_progress_tracker) %fresh start
end
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf)
parfor idx = 1:num_files
    %if mod(idx,1000) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %store parameters
    %file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{1};
    %just get all of them, baseline test. Everything that's not undecided
    valid_states = ~strcmpi(state_durations(:,end),'undecided');
    %state.count recorded in second col 
    state_durations = state_durations(valid_states,2);
    state_durations = cat(1,state_durations{:});
    %convert to time
    state_durations = state_durations * timestep;
    %ratecheck estimates
    Rcheck = curr_file.sim_results{4};
    %store durations, parameters, rate estimates
    file_data(idx,:) = {state_durations,curr_file.options,Rcheck};
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .1)) == 0 %at 10 percent
        progress = (progress / num_files) * 100;
        message = sprintf('----%.1f percent complete',progress);
        update_logfile(message,output_log)
    end
end
delete(gcp('nocreate'))

%search for jobs with identical parameters, collapse distributions
%get the randomized network parameters
job_params = cellfun(@(x)  [x.ItoE, x.EtoI],file_data(:,2),'UniformOutput',false);
job_params = vertcat(job_params{:});
uniq_params = unique(job_params,'rows');
num_jobs = size(uniq_params,1);
fprintf('----------------------\n')
fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)

%collapse duplicate job parameters
result_data = cell(num_jobs,3);
for idx = 1:num_jobs
    %find all matching
    curr_file = ismember(job_params,uniq_params(idx,:),'rows');
    %collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
    %rate estimates
    Rcheck_data = file_data(curr_file,3);
    curr_rate.Erate = mean(cellfun(@(x) x.Erate,Rcheck_data));
    curr_rate.Irate = mean(cellfun(@(x) x.Irate,Rcheck_data));
    result_data{idx,3} = curr_rate;
end


num_states = cellfun(@(x) numel(x),result_data(:,1));
mu_time = cellfun(@(x) mean(x),result_data(:,1));
fprintf('\n---jobs without data = %i\n',sum(num_states == 0))

Tmax = 300; %set a maximum duration for these plots 
Tmin = 1; %minimum duration 

overmax = mu_time > Tmax;
undermin = mu_time < Tmin;
Tinvalid = num_states == 0;
fprintf('\n---Tmax cutoff = %i\n',Tmax) %before log transform... 
fprintf('---%i / %i jobs above cutuff\n',sum(overmax),num_jobs)
fprintf('\n---Tmin cutoff = %i\n',Tmin)
fprintf('---%i / %i jobs under cutuff\n',sum(undermin),num_jobs)
fprintf('\n---%i / %i jobs without data\n',sum(Tinvalid),num_jobs)

Tinvalid = Tinvalid | overmax | undermin;

%---narrowing parameter range down
switch limit_prange
    case 'yes'
        %find parameters < value. x is otpions field ('EtoI') and y is cutoff value
        param_OOB = @(x,y) cellfun(@(z)  z.(x),result_data(:,2)) <= y;
        
        EI_valid = param_OOB('EtoI',EI_max);
        IE_valid = param_OOB('ItoE',IE_max);
        fprintf('\n---E-to-I cutoff = %.2f\n',EI_max)
        fprintf('---%i / %i jobs under cutuff\n',sum(EI_valid),num_jobs)
        fprintf('\n---I-to-E cutoff = %.2f\n',IE_max)
        fprintf('---%i / %i jobs under cutuff\n',sum(IE_valid),num_jobs)
        
        Tinvalid = Tinvalid | ~EI_valid | ~IE_valid;
end

fprintf('\nnum valid jobs = %i\n',sum(~Tinvalid))



%stats
switch outcome_stat
    case 'mu'
        outcome = cellfun(@(x) mean(x),result_data(:,1));
        Zlabel = 'mean duration (s)';
        figdir = fullfile(figdir,'mean_duration');
    case 'med'
        outcome = cellfun(@(x) median(x),result_data(:,1));
        Zlabel = 'median duration (s)';
        figdir = fullfile(figdir,'med_duration');
    case 'logmu'
        outcome = cellfun(@(x) mean(log10(x+eps)),result_data(:,1));
        Zlabel = {'mean stay duration';'log_{10}(s)'};
        figdir = fullfile(figdir,'logmean_duration');
        rad2keep = .15;
    case 'E-rate'
        outcome = cellfun(@(x) x.Erate,result_data(:,3));
        Zlabel = {'E spikerates';'Hz'};
        figdir = fullfile(figdir,'rates_excit');
        rad2keep = .15;
    case 'I-rate'
        outcome = cellfun(@(x) x.Irate,result_data(:,3));
        Zlabel = {'I spikerates';'Hz'};
        figdir = fullfile(figdir,'rates_inhib');
        rad2keep = .5;
end

switch limit_prange
    case 'yes'
        figdir = fullfile(figdir,'trimmed_range');
end

if ~isdir(figdir),mkdir(figdir);end


outcome = outcome(~Tinvalid);
%get the network parameters
EtoE = cellfun(@(x)  x.EtoE,result_data(~Tinvalid,2));
ItoE = cellfun(@(x)  x.ItoE,result_data(~Tinvalid,2));
EtoI = cellfun(@(x)  x.EtoI,result_data(~Tinvalid,2));
num_jobs = numel(outcome);


%save a big results file 
data_sz = whos('result_data');
data_sz = data_sz.bytes / 1e6; %in MB
if data_sz > 300
    fprintf('\nWARNING: summary file size = %.0f MB\n... summary file will not be saved\n',data_sz)
else
    fprintf('\nsaving summary file...')
    save(fullfile(resdir,'summary_file'),'result_data')
    fprintf('complete\n')
end

%surface plot
figure
fontsz = 20;
pointsz = 10; %was 30
plt_alph = .75;
num_gridlines = 400; 
xlin = linspace(min(ItoE),max(ItoE),num_gridlines);
ylin = linspace(min(EtoI),max(EtoI),num_gridlines);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(ItoE,EtoI,outcome);
f.Method = 'natural';
Z = f(X,Y);
%---if you need to smooth the plane---
filtSD = 2.5;
Z = imgaussfilt(Z,filtSD);
%---old mask implementation 
%mask = tril(ones(size(Z)),50);%mask = logical(flipud(mask));%Z(~mask) = NaN;
%rad2keep = .15;
mesh2keep = zeros(size(X));
for idx = 1:num_jobs
    x = ItoE(idx);y = EtoI(idx);z = outcome(idx);
    sphere_voxels = ((Y - y).^2 + (X - x).^2 + (Z - z).^2)<= rad2keep.^2;
    mesh2keep(sphere_voxels) = 1;
end
X(~mesh2keep) = NaN;Y(~mesh2keep) = NaN;Z(~mesh2keep) = NaN;
mesh(X,Y,Z,'FaceAlpha',plt_alph,'EdgeAlpha',plt_alph) 
axis tight; hold on
%these feel like they need to be slightly bigger 
%scatter3(ItoE,EtoI,outcome,25,'black','filled','MarkerEdgeAlpha',1,'MarkerEdgeAlpha',1) %give them outlines
scatter3(ItoE,EtoI,outcome,pointsz,'red','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1) 
xlabel({'within-pool inhibition';'(I-to-E strength)'},'FontWeight','b')
ylabel({'cross-pool inhibition';'(E-to-I strength)'},'FontWeight','b')
zlabel(Zlabel,'FontWeight','b')
set(gca,'FontSize',fontsz)
view(-45,27)
%view(-124,31)
%view(0,90) %that looks good too
hidden off
rotate3d on
%also see--
%https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html#bsovi2t
% switch rescale_plane
%     case 'on'
%         Zlim = get(gca,'ZLim');
%         Zlim(1) = 0;
%         Zlim(2) = max(outcome);
%         set(gca,'ZLim',Zlim);
%         caxis(Zlim)
% end
%savefig(fullfile(figdir,'surface_plot'))
%make it big
set(gcf,'units','normalized','outerposition',[0 0 .75 1])
%fix the labels 
xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
xh_pos = get(xh, 'Position');
%set(xh, 'Position',xh_pos+[-.05,.1,0],'Rotation',-24.5)
set(xh, 'Position',xh_pos+[.05,.1,0],'Rotation',20.5)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
yh_pos = get(yh, 'Position');
%set(yh, 'Position',yh_pos+ [.1,.15,0],'Rotation',13)
set(yh, 'Position',yh_pos+ [-.075,.15,0],'Rotation',-19)
set(gcf,'Renderer','painters')
print(fullfile(figdir,'surface_plot'),'-djpeg','-r400')%print high-res

%https://www.mathworks.com/matlabcentral/answers/41800-remove-sidewalls-from-surface-plots
savefig(gcf,fullfile(figdir,'surface_plot'),'compact')
%view(-24,24)

figure
hold off

Z = flipud(Z); %so the bottom left corner is low x, low y
% switch rescale_plane
%     case 'on'
%         Z(Z < 0) = 0;
%         Z(Z > max(outcome)) = max(outcome);
% end
%mask = tril(ones(size(Z)),100);
%mask = logical(mask);
%Z(~mask) = NaN;
colormap(parula)
imagesc(Z,'AlphaData',~isnan(Z))
xlabel('I-to-E strength')
ylabel('E-to-I strength')
title('State durations')
set(gca,'Fontsize',fontsz-4)
colb = colorbar;
colb.Label.String = Zlabel(end);
colb.FontSize = 16;
origXticks = get(gca,'Xtick');
Xticks = linspace(min(ItoE),max(ItoE),max(origXticks));
Xticks = Xticks(origXticks);
Xlabs = cellfun(@(x) sprintf('%.1f',x),num2cell(Xticks),'UniformOutput', false);
set(gca, 'XTickLabel', Xlabs)
origYticks = get(gca,'Ytick');
Yticks = linspace(max(EtoI),min(EtoI),max(origYticks));
%reset these real quick
origYticks = [1,origYticks(1:end-1)];
set(gca,'Ytick',origYticks);
Yticks = Yticks(origYticks);
Ylabs = cellfun(@(x) sprintf('%.2f',x),num2cell(Yticks),'UniformOutput', false);
set(gca, 'YTickLabel', Ylabs')
set(gca,'FontSize',fontsz-4)
print(fullfile(figdir,'heatmap'),'-djpeg')
savefig(gcf,fullfile(figdir,'heatmap'),'compact')



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


hold off; close all

% 
% %try countour plot
% xlin = linspace(min(ItoE),max(ItoE),num_jobs);
% ylin = linspace(min(EtoI),max(EtoI),num_jobs);
% [X,Y] = meshgrid(xlin,ylin);
% f = scatteredInterpolant(ItoE,EtoI,outcome);
% f.Method = 'natural';
% Z = f(X,Y);
% switch rescale_plane
%     case 'on'
%         Z(Z < 0) = 0;
%         Z(Z > max(outcome)) = max(outcome);
% end
% mask = tril(ones(size(Z)),100);
% mask = logical(flipud(mask));
% Z(~mask) = NaN;
% figure
% [c,h] = contour(X,Y,Z,'LineWidth',2);
% levs = h.LevelList;
% h.LevelList = levs(1:2:end);
% clabel(c,h);
% title('State durations')
% xlabel('I-to-E strength')
% ylabel('E-to-I strength')
% set(gca,'FontSize',fontsz)
% print(fullfile(figdir,'contour'),'-djpeg')
% 
% 
% %this 3D contour plot is really good actually
% H = figure;
% [c,h] = contour3(X,Y,Z,'LineWidth',5);
% levs = h.LevelList;
% h.LevelList = levs(1:2:end);
% %clabel(c,h);
% axis tight; hold on
% plot3(ItoE,EtoI,outcome,'.','MarkerSize',25,'color','red') %nonuniform
% xlabel('I-to-E strength')
% ylabel('E-to-I strength')
% zlabel(Zlabel)
% set(gca,'FontSize',fontsz)
% savefig(H,fullfile(figdir,'contour3D_plot'),'compact')
% 
% 
% hold off
% figure
% histogram(outcome)
% ylabel('frequency')
% xlabel(Zlabel)
% title('State durations')
% set(gca,'FontSize',fontsz)
% print(fullfile(figdir,'histogram'),'-djpeg')

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
