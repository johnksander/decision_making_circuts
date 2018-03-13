clear
clc
format compact
hold off;close all

rescale_plane = 'on';
outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'

%result summaries
fontsz = 16;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_stims';
basedir = '/home/acclab/Desktop/ksander/rotation';
figdir = fullfile(basedir,'Results',[sim_name '_figures']);
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_BL'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
network_pair_info = load(fullfile(basedir,'helper_functions','network_pairs'));
network_pair_info = network_pair_info.network_pairs;
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],...
    network_pair_info,'UniformOutput',false);
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],...
    network_pair_info,'UniformOutput',false);

%get results
file_data = cell(num_files,2);
for idx = 1:num_files
    curr_file = load(output_fns{idx});
    %store parameters
    file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{:};
    %take only stimulus state durations
    state_durations = state_durations{1}(:,1);
    %     state_durations = cellfun(@(x) x(:,1),state_durations,'UniformOutput',false);
    %     state_durations = vertcat(state_durations{:});
    state_durations = vertcat(state_durations{:}); %ooo that's annoying
    %convert to time
    state_durations = state_durations * timestep;
    file_data{idx,1} = state_durations;
end

%search for jobs with identical parameters, collapse distributions
%get the randomized network parameters
job_params = cellfun(@(x)...
    [x.ItoE, x.EtoI,unique(x.trial_stimuli),find(strcmpi(x.stim_targs, stimtarg_vals))],...
    file_data(:,2),'UniformOutput',false); %matching "network_pair_info" format
job_params = vertcat(job_params{:});
uniq_params = unique(job_params,'rows');
num_jobs = size(uniq_params,1);
fprintf('----------------------\n')
fprintf('num jobs = %i\nunique parameter sets = %i\nduplicates = %i\n',num_files,num_jobs,num_files - num_jobs)

%collapse duplicate job parameters
result_data = cell(num_jobs,2);
Nruns = NaN(num_jobs,1); %record the number of successful jobs..
for idx = 1:num_jobs
    %find all matching
    curr_file = ismember(job_params,uniq_params(idx,:),'rows');
    Nruns(idx) = sum(curr_file);
    fprintf('---\nparameter set %.3f %.3f %.3f %s, n files = %i\n',uniq_params(idx,1:3),...
        stimtarg_vals{uniq_params(idx,4)},Nruns(idx))
    %collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    fprintf('------n states = %i\n',numel(result_data{idx,1}))
    %just grab the first options file... that shouldn't matter here
    result_data{idx,2} = file_data{find(curr_file,1),2};
end

need_more = cellfun(@(x) numel(x{1}),num2cell(result_data,2));
need_more = need_more < 10000;

%just grab some simple stats here
num_states = cellfun(@(x) numel(x),result_data(:,1));
mu_duration = cellfun(@(x) mean(x),result_data(:,1));
med_duration = cellfun(@(x) median(x),result_data(:,1));
logmu_dur = cellfun(@(x) mean(log(x)),result_data(:,1));
%get the network parameters
EtoE = cellfun(@(x)  x.EtoE,result_data(:,2));
ItoE = cellfun(@(x)  x.ItoE,result_data(:,2));
EtoI = cellfun(@(x)  x.EtoI,result_data(:,2));
stim_value = cellfun(@(x)  unique(x.trial_stimuli),result_data(:,2));

net_type = num2cell(num2cell(uniq_params),2);
net_type = cellfun(@(x)  [x(1:3),stimtarg_vals{x{4}}],net_type,'UniformOutput',false);
net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],net_type,'UniformOutput',false);
net_type = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],net_type,'UniformOutput',false);

switch outcome_stat
    case 'mu'
        outcome = mu_duration;
        %Zlabel = 'mean duration (s)';
        Zlabel = 'state duration (s)';
        figdir = fullfile(figdir,'mean_duration');
    case 'med'
        outcome = med_duration;
        Zlabel = 'median duration (s)';
        figdir = fullfile(figdir,'med_duration');
    case 'logmu'
        outcome = logmu_dur;
        %Zlabel = 'mean log duration [ log(s) ]';
        Zlabel = 'log(s) state duration';
        figdir = fullfile(figdir,'logmean_duration');
end

if ~isdir(figdir),mkdir(figdir);end



num_types = numel(network_pair_info);
plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        hold on
        %find the right results for network set-up
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = cell2mat(result_data(curr_data,1));
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = log(curr_data);
                    [kde,kde_i] = ksdensity(curr_data);
                    plot(kde_i,kde,'LineWidth',2);
                case 'mu'
                    histogram(curr_data)
            end
        end
        %take control data
        BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
        curr_data = cellfun(@(x) isequal(x,BLinfo),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = cell2mat(result_data(curr_data,1));
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = log(curr_data);
                    [kde,kde_i] = ksdensity(curr_data);
                    plot(kde_i,kde,'LineWidth',2);
                case 'mu'
                    histogram(curr_data)
            end
        end
        hold off
        legend(sprintf('%.0fHz',curr_net_info{j,3}))
        
        if plt_idx == 9 || plt_idx == 10
            xlabel(Zlabel)
        end
        if plt_idx == 1 || plt_idx == 2
            title(sprintf('%s networks',curr_net_info{j,4}),'Fontsize',14)
        end
        if mod(plt_idx,2) == 1
            ylabel(sprintf('network #%i p(x)',idx))
        end
        
    end
end
orient tall
linkaxes(h,'x')
axis tight
print(fullfile(figdir,'duration_dists'),'-djpeg')



num_states = cellfun(@(x) numel(x{1}),num2cell(result_data,2));
need_more = num_states < 10000;
need_more = net_type(need_more);
current_count = num_states(num_states < 10000);
for idx = 1:num_types
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),need_more,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        if sum(curr_data) > 0
           fprintf('\nnetwork #%i %s has < 10k states (%i)',idx,curr_net_info{j,4},current_count(curr_data))
           fprintf('\n---parameter set:\n')
           disp(curr_net_info(j,:))
           fprintf('\n------------------------\n')
        end
        
        BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
        curr_data = cellfun(@(x) isequal(x,BLinfo),need_more,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        if sum(curr_data) > 0
           fprintf('\nnetwork #%i %s has < 10k states (%i)',idx,BLinfo{4},current_count(curr_data))
           fprintf('\n---parameter set:\n')
           disp(BLinfo)
           fprintf('\n------------------------\n')
        end
        
        
    end
        
end

maxout = (150*10)+3000;
%stims
lol = 3002:10:maxout;
%BL
lol = 3001:10:maxout;
lol = 3003:10:maxout;
lol = 3005:10:maxout;
lol = 3007:10:maxout;
lol = 3009:10:maxout;

%hacky code for irregular jobs

%BLname = 'parsweep_stims_BL';
%BL_fns = dir(fullfile(basedir,'Results',BLname,['*',BLname,'*.mat']));

% fprintf('\nWARN: hardcode-resetting parameters\n')
% for idx = 1:numel(network_pair_info)
%    fix_these = strcmpi(network_pair_info{idx},'fast');
%    fix_these = [fix_these(:,2:end),[0;0]];
%    network_pair_info{idx}(logical(fix_these)) = {200};
% end





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
keyboard





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
print(fullfile(figdir,'net_inhib_scatter'),'-djpeg')
keyboard


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
figure;
mesh(X,Y,Z) %interpolated
axis tight; hold on
plot3(ItoE,EtoI,outcome,'.','MarkerSize',25,'color','red') %nonuniform
xlabel('I-to-E strength')
ylabel('E-to-I strength')
zlabel(Zlabel)
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
set(gca,'FontSize',fontsz)
savefig(fullfile(figdir,'surface_plot'))



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
figure
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
