clear
clc
format compact
hold off;close all
%note: if you ever wana make this nice... make a single function for
%creating a paramater ID cell array (or table) from an array of options structure.
%this function would be used for the template nets from get_network_params() and result files.
%Put this same function in network_spiking_results, etc. Also use cellfun(@(x) isequal(x,table2cell(curr_net_info(j,:))),Psets)

opt = struct();
opt.print_anything = 'yes'; %'yes' | 'no';
opt.multiple_stimuli = 'no';
opt.valid_states = 'stay'; %'stay' | 'all'; undecided is always invalid, 'all' gives stay & leave
opt.outcome_stat = 'mu';  %'mu' | 'med' | 'logmu'
opt.pulse_stim = 'off'; %'yes' | 'total_time' | 'rem' | 'off' whether to treat durations as samples (rem = time during sample)
opt.parfor_load = 'on'; %on/off, must also (un)comment the actual for... line
opt.params2match = {'conn','stim'}; %!!!IMPORTANT!!! specify how results are matched to network types
%this can be at most {'conn','stim'}. That specifies matching on connection strengths, stimulus values


%specify simulations
Snames = {'testing_dt02','testing_dt10'};
figdir = cellfun(@(x) sprintf('figures_%s',x),Snames,'UniformOutput',false);

basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

Njobs = numel(Snames);
data = cell(Njobs,1);
job_params = cell(Njobs,1);
for idx = 1:numel(Snames)
    [data{idx},job_params{idx}] = get_data(basedir,Snames{idx},figdir{idx},opt);
end

sample_sz = 5e3; %seconds
dt_job = cellfun(@(x) x.timestep,job_params);
Nfiles = cellfun(@numel,data); %files per job 
Tjob = cellfun(@(x) x.tmax,job_params); %tmax 
Ndraw = ceil(sample_sz ./ Tjob); %how many files to draw
Tsamp = Tjob .* Ndraw; %how long is each sample, total (s)
Nsamp = 10e3; %how many total samplings 
results = NaN(Nsamp,2);

samp_inds = @(n,k) randi(n,k,1); %without replacement 
cell2vec = @(x) cat(1,x{:});

num_workers = 24;
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_stay_durations')})
special_progress_tracker = fullfile(basedir,'SPT.txt');
if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start

for idx = 1:Nsamp
    
    %sample indicies 
    inds = arrayfun(samp_inds,Nfiles,Ndraw,'UniformOutput',false);
    sample = cellfun(@(x,i) x(i),data,inds,'UniformOutput',false);
    sample = cellfun(cell2vec,sample,'UniformOutput',false);    
    sample = cellfun(@mean,sample);
    results(idx,:) = sample;
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(Nsamp * .1)) == 0 %at 10 percent
        progress = (progress / Nsamp) * 100;
        fprintf('----%.1f percent complete\n',progress);
    end
end

delete(gcp('nocreate'))
delete(special_progress_tracker)

results = num2cell(results,1);

labs = arrayfun(@(x) sprintf('\\Deltat = %.2f ms',x*1e3),dt_job,'Uniformoutput',false);

figure;orient tall
subplot(2,1,1)
scatter(results{1},results{2})
%axis tight
lims = axis;
lims(1:2:end) = min(lims);
lims(2:2:end) = max(lims);
axis(lims)
title(sprintf('mean stay duration (%is sample)',unique(Tsamp)))
xlabel(labs{1});ylabel(labs{2})
xtickformat('%is')
ytickformat('%is')
%show correlation
rho = corr(results{1},results{2},'type','spearman');
set(gca,'FontSize',12)
text(.05,.1,sprintf('\\rho = %.4f',rho),'Units','normalized','FontWeight','bold','FontSize',13)


subplot(2,1,2)
histogram(results{1}-results{2})
xlabel(sprintf('sample mean difference\n[\\Deltat  %.2f ms] - [\\Deltat  %.2f ms]',...
    dt_job(1)*1e3,dt_job(2)*1e3))
ylabel('frequency')
xtickformat('%is')
set(gca,'FontSize',12)

FN = sprintf('dt%i_dt%i_sample_comparisons',dt_job(1)*1e5,dt_job(2)*1e5);
print(FN,'-djpeg')
close all

%options.timestep = .1e-3; 
% 
% R = 2e3; %Hz
% L = R.*dt_job;
% T = 1;
% N = 200;
% Nt = T./ dt_job;
% test = cell(Njobs,1);
% for idx = 1:Njobs
%     curr_Nt = Nt(idx);
%     spikes = zeros(N,1);
%     for tidx = 1:curr_Nt
%         spikes = spikes + poissrnd(L(idx),N,1);
%     end
%     test{idx} = spikes;
% end
% figure
% histogram(test{1});hold on
% histogram(test{2});hold off
% legend(labs,'Location','best','box','off')

return
%specify simulations
Snames = {'testing_dt02','testing_dt10','testing_dt25'};
figdir = cellfun(@(x) sprintf('figures_%s',x),Snames,'UniformOutput',false);

Njobs = numel(Snames);
data = cell(Njobs,1);
job_params = cell(Njobs,1);
for idx = 1:numel(Snames)
    [data{idx},job_params{idx}] = get_data(basedir,Snames{idx},figdir{idx},opt);
end

data = cellfun(cell2vec,data,'UniformOutput',false);
dt_job = cellfun(@(x) x.timestep,job_params);


figure; orient tall
subplot(2,1,1)
plot(dt_job*1e3,cellfun(@mean,data),'-o','LineWidth',2)
xlim([0,max(get(gca,'Xlim')) + .02])
xtickformat('%.2f')
xlabel('\Deltat (ms)')
ylabel('duration (s)')
set(gca,'FontSize',14)

subplot(2,1,2)
plot(dt_job*1e3,cellfun(@(x) mean(log10(x)),data),'-o','LineWidth',2)
xlim([0,max(get(gca,'Xlim')) + .02])
xtickformat('%.2f')
xlabel('\Deltat (ms)')
ylabel('duration log_{10}(s)')
set(gca,'FontSize',14)
print('dt_duration_lineplot','-djpeg')



return

%equate X axes for all figs of the same type
%Snames = Snames(1:end-1);
figdir = 'figures_timestep_testing';
figdir = fullfile(basedir,'Results',figdir,'durations');
savedir = fullfile(figdir,'Xsynced');
if ~isdir(savedir),mkdir(savedir);end
ftypes = {'total_time','total_time_log'};
%ftypes = {'decision_timing_log','decision_timing','total_time','total_time_log','total_samples'};

invid_cols = {'no','no'}; %equate x axes within column (e.g. slow vs fast)

for idx = 1:numel(ftypes)
    
    Xls = NaN(numel(Snames),2);
    switch invid_cols{idx}
        case 'yes'
            Xls = cat(3,Xls,Xls);
    end
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        switch invid_cols{idx}
            case 'yes'
                ax = num2cell(findall(h, 'type', 'axes'));
                %left column even, right column odd
                ax = cellfun(@(x) x.XLim,ax,'UniformOutput',false);
                ax = [ax(2:2:end),ax(1:2:end)]; %sort to columns
                ax = num2cell(ax,1);
                ax = cellfun(@cell2mat,ax,'UniformOutput',false);
                ax = cellfun(@(x) [min(x(:,1)),max(x(:,2))],ax,'UniformOutput',false);
                Xls(fidx,:,:) = cat(3,ax{:}); %put in 3rd D
            otherwise
                Xls(fidx,:) = h.CurrentAxes.XLim;
        end
        close all
    end
    switch invid_cols{idx}
        case 'yes'
            xlims = [min(Xls(:,1,:)),max(Xls(:,2,:))];
            xlims = squeeze(xlims); %now Xlims column is Xlims for a column
        otherwise
            xlims = [min(Xls(:,1)),max(Xls(:,2))];
    end
    %xlims = [0,.8];
    for fidx = 1:numel(Snames)
        fn = fullfile(figdir,[Snames{fidx} '_' ftypes{idx} '.fig']);
        h = openfig(fn);
        ax = findobj(h,'Type','Axes');
        switch invid_cols{idx}
            case 'yes' %left column even, right column odd
                set(ax(2:2:end),'XLim',xlims(:,1))
                set(ax(1:2:end),'XLim',xlims(:,2))
            otherwise
                set(ax,'XLim',xlims)
        end
        ch = allchild(ax);
        if iscell(ch)
            ch = cat(1,ch{:});
        end
        %set(ch,'Normalization','probability')
        %binwit = cellfun(@(x) x.XLim,num2cell(ax),'UniformOutput',false);
        %binwit = cellfun(@range,binwit) ./ 50;
        %binwit = num2cell(binwit);
        %nodata = cellfun(@(x) isempty(x.Children),num2cell(ax));
        %cellfun(@(x,y) set(x,'BinWidth',y),num2cell(ch),binwit(~nodata));   %set(ch,'BinWidth',.1)
        Fdir = fullfile(savedir,ftypes{idx});
        leg = get(gca,'Legend');
        set(leg,'location','best')
        if ~isdir(Fdir),mkdir(Fdir);end
        print(fullfile(Fdir,[Snames{fidx} '_' ftypes{idx}]),'-djpeg')
        close all
    end
    
end



function [file_data,job_params] = get_data(basedir,sim_name,figdir,opt)
hold off;close all

outcome_stat = opt.outcome_stat;
pulse_stim = opt.pulse_stim;
print_anything = opt.print_anything; %'yes' | 'no';
summary_stats = 'no'; %summary_stats = opt.summary_stats;
params2match = opt.params2match;
figdir = fullfile(basedir,'Results',figdir,'durations');
resdir = fullfile(basedir,'Results',sim_name);
%output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = dir(fullfile(resdir,'*.mat')); %use this for unrestricted loading
warning('loading all mat files from results directory!!')
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
%checking for previously saved data
svdir = fullfile(figdir,'data');if ~isdir(svdir),mkdir(svdir);end
svFN = [sim_name '_%s.mat'];
switch pulse_stim
    case 'yes'
        svFN = sprintf(svFN,'total_samples');
    case 'rem'
        svFN = sprintf(svFN,'decision_timing');
    otherwise
        svFN = sprintf(svFN,'total_time');
end
if exist(fullfile(svdir,svFN)) > 0,load_summary = true;else,load_summary = false;end

switch opt.multiple_stimuli
    case 'yes'
        param_varnams = {'ItoE','EtoI','stim_A','stim_B','targ_cells'};
        %line below may be important
        %IDvars = param_varnams(~ismember(param_varnams,'stim_B')); %stim B not particular to network type
    case 'no'
        param_varnams = {'ItoE','EtoI','stim','targ_cells'};
end
%for indexing the result paramters
IDvars = [];
if sum(strcmp('conn',params2match)) > 0,IDvars = {'ItoE','EtoI'};end
if sum(strcmp('stim',params2match)) > 0,IDvars = [IDvars,param_varnams(startsWith(param_varnams,'stim'))];end


%get general options file from the first file
gen_options = load(output_fns{1});
gen_options = gen_options.options;
timestep = gen_options.timestep;

switch pulse_stim
    case 'off'
        %skip this business
    otherwise
        %pulse duration... kinda hardcoded here
        error('get this from  options dude, was previously striped from sim name')
end

num_files = numel(output_fns);
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'};

%info on the specific network parameters in this simulation
% num_net_types = 10;
% num_pairs = 5;
% pair_inds = num2cell(reshape(1:num_net_types,[],num_pairs)); %just gives a cell array for pair indicies
% network_pair_info = cell(num_pairs,1);
% for idx = 1:num_pairs
%     curr_params = cellfun(@(x) get_network_params(x,gen_options),pair_inds(:,idx),'UniformOutput',false);
%     switch opt.multiple_stimuli
%         case 'yes'
%             curr_params = cellfun(@(x)...
%                 [x.ItoE, x.EtoI,num2cell(x.trial_stimuli),x.stim_targs],...
%                 curr_params,'UniformOutput',false); %matching "network_pair_info" format
%         otherwise
%             curr_params = cellfun(@(x)...
%                 {x.ItoE, x.EtoI,unique(x.trial_stimuli), x.stim_targs},...
%                 curr_params,'UniformOutput',false); %matching "network_pair_info" format
%     end
%     curr_params = cat(1,curr_params{:});
%     T = cell2table(curr_params,'VariableNames',param_varnams);
%     curr_types = T.targ_cells;
%     curr_types = strrep(curr_types,'Estay','fast'); curr_types = strrep(curr_types,'Eswitch','slow');
%     T.Properties.RowNames = curr_types;
%     network_pair_info{idx} = T;
% end

fprintf('\n---loading simulation: %s\n',sim_name)

if ~load_summary
    %get results
    switch opt.parfor_load
        case 'off'
            fprintf('\nparfor disabled\n')
        case 'on'
            num_workers = 24;
            c = parcluster('local');
            c.NumWorkers = num_workers;
            parpool(c,c.NumWorkers,'IdleTimeout',Inf,'AttachedFiles',{which('find_stay_durations')})
            special_progress_tracker = fullfile(basedir,'SPT.txt');
            if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start
    end
    switch opt.valid_states %select states for analysis
        case 'all'
            warning('\nfind_stay_durations() disabled, all non-undecided states')
    end
    
    file_data = cell(num_files,2);
    parfor idx = 1:num_files
        %for idx = 1:num_files
        switch opt.parfor_load
            case 'off'
                if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
        end
        curr_file = load(output_fns{idx});
        %get state durations
        state_durations = curr_file.sim_results;
        state_durations = state_durations{1};
        switch opt.valid_states %select states for analysis
            case 'all'
                %just get all of them, baseline test. Everything that's not undecided
                keep_states = ~strcmpi(state_durations(:,end),'undecided');
                %state.count recorded in second col
                state_durations = state_durations(keep_states,2);
                state_durations = cat(1,state_durations{:});
                %convert to time
                state_durations = state_durations * timestep;
            case 'stay'
                state_durations = find_stay_durations(state_durations,curr_file.options,'verify');
                switch pulse_stim
                    case 'yes' %just do this now while options is handy
                        state_durations = state_durations.samples;
                    case 'rem' %look at when IN the sample switch happened
                        state_durations = state_durations.decision_time;
                    case 'total_time'
                        state_durations = state_durations.duration;
                    case 'off'
                        state_durations = state_durations.duration;
                end
        end
        
        %store durations & parameters
        file_data(idx,:) = {state_durations,curr_file.options};
        
        switch opt.parfor_load
            case 'on'
                progress = worker_progress_tracker(special_progress_tracker);
                if mod(progress,floor(num_files * .05)) == 0 %at half a percent
                    progress = (progress / num_files) * 100;
                    fprintf('----%.1f percent complete\n',progress);
                end
        end
    end
    switch opt.parfor_load
        case 'on'
            delete(gcp('nocreate'))
            delete(special_progress_tracker)
    end
    
    fprintf('\nsaving data...\n')
    save(fullfile(svdir,svFN),'file_data')
elseif load_summary
    fprintf('\nloading saved summary data...\n')
    file_data = load(fullfile(svdir,svFN));
    file_data = file_data.file_data;
end

%take params from first options file 
job_params = file_data{1,2};

file_data = file_data(:,1);
%file_data = cat(1,file_data{:});
end


