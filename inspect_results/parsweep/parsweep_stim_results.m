clear
clc
format compact
hold off;close all

rescale_plane = 'on';
outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu' 
pulse_stim = 'rem'; %'yes' | 'no' | 'rem' whether to treat durations as samples (rem = time during sample)
print_anything = 'yes'; %'yes' | 'no';

%result summaries
fontsz = 16;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%specify simulation
%---sim setup-----------------
%run_num = 3; %for seperate figure directories & restricting file loading 
sim_name = 'network_spiking_P2_1';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))
figdir = fullfile(basedir,'Results','network_spiking_pulse_figures','durations',sim_name);
%figdir = fullfile(basedir,'Results',[sim_name '_figures']);
%figdir = fullfile(figdir,sprintf('run%i_figs',run_num));
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
%output_fns = dir(fullfile(resdir,sprintf('*%s*run%i*.mat',sim_name,run_num)));
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_BL'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
%this file was being created "by hand" with parsweep_find_examples
%network_pair_info = load(fullfile(basedir,'helper_functions','network_pairs'));
%network_pair_info = network_pair_info.network_pairs;
%take this stuff directly from get_network_params()
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

%this code is holdover from when it was loaded.. also dumb. 
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')],...
    network_pair_info,'UniformOutput',false);
network_pair_info = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')],...
    network_pair_info,'UniformOutput',false);

%get results
file_data = cell(num_files,2);
for idx = 1:num_files
    if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %store parameters
    file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{1};
    %take only stimulus state durations
    state_durations = state_durations{1}(:,1);
    %     state_durations = cellfun(@(x) x(:,1),state_durations,'UniformOutput',false);
    %     state_durations = vertcat(state_durations{:});
    state_durations = cat(1,state_durations{:}); %ooo that's annoying
    %convert to time
    state_durations = state_durations * timestep;
    switch pulse_stim
        case 'yes' %just do this now while options is handy 
            state_durations = floor(state_durations ./ sum(curr_file.options.stim_pulse));
        case 'rem' %look at when IN the sample switch happened
            state_durations = mod(state_durations,sum(curr_file.options.stim_pulse));
    end
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


switch pulse_stim
    case 'yes' 
        unit_measure = 'samples';
    case 'rem' 
         unit_measure = 's - onset';
    otherwise
        unit_measure = 's'; %like "seconds" not samples
end

switch outcome_stat
    case 'mu'
        outcome = mu_duration;
        %Zlabel = 'mean duration (s)';
        Zlabel = sprintf('state duration (%s)',unit_measure);
        figdir = fullfile(figdir,'mean_duration');
    case 'med'
        outcome = med_duration;
        Zlabel =  sprintf('median duration (%s)',unit_measure);
        figdir = fullfile(figdir,'med_duration');
    case 'logmu'
        outcome = logmu_dur;
        %Zlabel = 'mean log duration [ log(s) ]';
        Zlabel = sprintf('log(%s) state duration',unit_measure);
        figdir = fullfile(figdir,'logmean_duration');
end

switch pulse_stim
    case {'yes','rem'}
        figdir = fullfile(figdir,'measured_in_samples');
end


if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;
alph = .5;

num_types = numel(network_pair_info);
plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{idx};
    for j = 1:2
        
        plt_idx = plt_idx + 1;
        h(plt_idx) = subplot(5,2,plt_idx);
        hold on
        
        
        %get the right color
        if strcmpi(curr_net_info{j,end},'slow')
            lcol = matorange;
        elseif strcmpi(curr_net_info{j,end},'fast')
            lcol = matblue;
        end
        
        %find the right results for network set-up
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        curr_data = cell2mat(result_data(curr_data,1));
                
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case 'logmu'
                    curr_data = curr_data(curr_data ~= 0);%inf errors
                    curr_data = log(curr_data);
            end
            switch pulse_stim
                case 'yes'
                    histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
                otherwise
                    [kde,kde_i] = ksdensity(curr_data);
                    area(kde_i,kde,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
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
                    curr_data = curr_data(curr_data ~= 0);%inf errors
                    curr_data = log(curr_data);
            end
            switch pulse_stim
                case 'yes'
                    histogram(curr_data,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
                otherwise
                    [kde,kde_i] = ksdensity(curr_data);
                    area(kde_i,kde,'FaceColor',lcol,'EdgeColor',lcol,'FaceAlpha',alph);
            end  
        end
        hold off
        legend(sprintf('%.0fHz',curr_net_info{j,3}),'location','best')
        
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

fig_fn = 'duration_dists';

switch pulse_stim
    case 'rem' %rename for special case 
        fig_fn = 'duration_after_stim_onset';
end
switch print_anything
    case 'yes'
        print(fullfile(figdir,fig_fn),'-djpeg')
end

%info about simulation
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

%summary stats
fprintf('\n------------------------\n')
fprintf('Summary statistics\n')
fprintf('%s:\n',Zlabel)
fprintf('------------------------\n\n')
%outcome = logmu_dur;
keyboard
for idx = 1:num_types
    fprintf('\n------------------------\nNetwork #%i\n',idx)
    curr_net_info = network_pair_info{idx};
    Xstim = cell(2,1); %for testing difference between stim distributions
    for j = 1:2
        fprintf('---type: %s\n',curr_net_info{j,end})
        %find the right results for network set-up
        curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        Xstim{j} = result_data{curr_data,1};
        fprintf('            stimulus = %.2f\n',outcome(curr_data))
        
        %take control data
        BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
        curr_data = cellfun(@(x) isequal(x,BLinfo),net_type,'UniformOutput',false);
        curr_data = cat(1,curr_data{:});
        fprintf('            baseline = %.2f\n',outcome(curr_data))  
    end
    
    fprintf('\n---hyp. test: mu stim durrations\n')
    switch outcome_stat
        case 'logmu'
            Xstim = cellfun(@log,Xstim,'UniformOutput',false);
    end
    [~,pval] = ttest2(Xstim{1},Xstim{2}); %regular old t-test
    fprintf('t-test p = %.3f\n',pval)
    [CI,H] = boot_mean_diffs(Xstim{1},Xstim{2},10000);
    fprintf('bootstrap test: %s\n',H)
    fprintf('bootstrap CI: %.2f, %.2f\n',CI)
    
    
end
keyboard

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
switch print_anything
    case 'yes'
        print(fullfile(figdir,'net_inhib_scatter'),'-djpeg')
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
switch print_anything
    case 'yes'
        print(fullfile(figdir,'heatmap'),'-djpeg')
end

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
switch print_anything
    case 'yes'
        print(fullfile(figdir,'contour'),'-djpeg')
end

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
switch print_anything
    case 'yes'
        savefig(H,fullfile(figdir,'contour3D_plot'),'compact')
end

hold off
figure
histogram(outcome)
ylabel('frequency')
xlabel(Zlabel)
title('State durations')
set(gca,'FontSize',fontsz)
switch print_anything
    case 'yes'
        print(fullfile(figdir,'histogram'),'-djpeg')
end


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
