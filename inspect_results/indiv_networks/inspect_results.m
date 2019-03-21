clear
clc
format compact
hold off;close all

outcome_stat = 'logmu';  %'mu' | 'med' | 'logmu' ||| 'E-rate' | 'I-rate'

%specify simulation
%---sim setup-----------------
sim_name = 'test_real_durations';
jobID = 3; %set to NaN for all jobs
basedir = '/home/acclab/Desktop/ksander/rotation/project';
figdir = fullfile(basedir,'Results',['figures_' sim_name]);
resdir = fullfile(basedir,'Results',sim_name);
addpath(fullfile(basedir,'helper_functions'))
%result summaries
fontsz = 12;
stim_labels = {'stim A','stim B'};

%get results
if ~isnan(jobID)
    figdir = fullfile(figdir,sprintf('jobID_%i',jobID)); %make job specific directory 
    output_fns = sprintf('%s_%i.mat',sim_name,jobID);
    output_fns = dir(fullfile(resdir,['*',output_fns]));
else
    output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
end
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
num_files = numel(output_fns);
file_data = cell(num_files,3);

%get the options struct from first file 
options = load(output_fns{1});
options = options.options;
timestep = options.timestep;

for idx = 1:num_files
    %if mod(idx,1000) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %store parameters
    %file_data{idx,2} = curr_file.options;
    %get state durations
    state_durations = curr_file.sim_results;
    state_durations = state_durations{1};
    %just get all stay states. Everything that's not undecided
    valid_states = startsWith(state_durations(:,end),'stim');
    %state.count recorded in second col 
    state_durations = state_durations(valid_states,2);
    state_durations = cat(1,state_durations{:});
    %convert to time
    state_durations = state_durations * timestep;
    file_data{idx,1} = state_durations;
    %     %ratecheck estimates
    %     Rcheck = curr_file.sim_results{4};
    %     %store durations, parameters, rate estimates
    %     file_data(idx,:) = {state_durations,curr_file.options,Rcheck};
    
end

state_durations = file_data(:,1);
state_durations = cat(1,state_durations{:});

%get the options struct from first file  
if ~isdir(figdir),mkdir(figdir);end


num_states = numel(state_durations);
stim = unique(options.trial_stimuli);

figure
subplot(2,1,1)
histogram(state_durations,num_states)
title(sprintf('Network:    E-I = %.2f,    I-E = %.2f\nstim = %.1fHz, N states = %i',options.EtoI,options.ItoE,stim,num_states))
xlabel('duration (s)');ylabel('frequecy');legend(sprintf('\\mu = %.2f',mean(state_durations)),'Location','best')
set(gca,'FontSize',fontsz)
subplot(2,1,2)
histogram(log10(state_durations),num_states)
xlabel('duration log_{10}(s)');ylabel('frequecy');legend(sprintf('\\mu = %.2f',mean(log10(state_durations))),'Location','best')
set(gca,'FontSize',fontsz)
print(fullfile(figdir,'state_durations'),'-djpeg')
