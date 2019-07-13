clear
clc
format compact
hold off;close all


%This code is intended to move result files that do not match some criteria
%(e.g. state duration time, parameter values, etc) into a subdirectory.
%This preserves all results from a parameter sweep, while simplifying
%analyses after determining those criteria. 

sim_name = 'parsweep_D2t_baseline'; %MUST CHANGE CRITERIA TOO

%----CRITERIA----
%connection strengths
%EI_max = .75; IE_max = 8; %these set connection maximums
%for the fastD
EI_max = .75; IE_max = 12.5; %these set connection maximums

%state durations
Tmax = 300; %set a maximum duration for these plots 
Tmin = 1; %minimum duration 
min_states = 1; %must have at least 1 state duration... 



num_workers = 24; %for parfor loading results 
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..

%specify simulation
%---sim setup-----------------
basedir = '/home/acclab/Desktop/ksander/rotation/project';
resdir = fullfile(basedir,'Results',sim_name);
inval_dir = fullfile(resdir,'invalid_results');
addpath(fullfile(basedir,'helper_functions'))
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = {output_fns.name};
if ~isdir(inval_dir),mkdir(inval_dir);end
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
    curr_file = load(fullfile(resdir,output_fns{idx}));
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
delete(gcp('nocreate'));delete(special_progress_tracker);delete(output_log)

%this code was adapted from the code that collapses identical job
%parameters. Not efficient this way, but I'm not spending the time recoding stuf... 

result_data = cell(num_files,3);
for idx = 1:num_files
    %%find all matching
    %curr_file = ismember(job_params,uniq_params(idx,:),'rows');
    curr_file = idx;
    %was: collapse & reallocate
    result_data{idx,1} = cell2mat(file_data(curr_file,1));
    %just grab the options struct
    result_data{idx,2} = file_data{curr_file,2};
    %rate estimates
    Rcheck_data = file_data(curr_file,3);
    curr_rate.Erate = mean(cellfun(@(x) x.Erate,Rcheck_data));
    curr_rate.Irate = mean(cellfun(@(x) x.Irate,Rcheck_data));
    result_data{idx,3} = curr_rate;
end


num_states = cellfun(@(x) numel(x),result_data(:,1));
mu_time = cellfun(@(x) mean(x),result_data(:,1));
fprintf('\n---jobs without data = %i\n',sum(num_states == 0))

overmax = mu_time > Tmax;
undermin = mu_time < Tmin;
underN = num_states < min_states;
fprintf('\n---Tmax cutoff = %i\n',Tmax) %before log transform... 
fprintf('---%i / %i jobs above cutoff\n',sum(overmax),num_files)
fprintf('\n---Tmin cutoff = %i\n',Tmin)
fprintf('---%i / %i jobs under cutoff\n',sum(undermin),num_files)
fprintf('\n---N states cutoff < %i\n',min_states)
fprintf('---%i / %i jobs under cutoff\n',sum(underN),num_files)

Tinvalid = underN | overmax | undermin;

%---narrowing parameter range down

%find parameters < value. x is otpions field ('EtoI') and y is cutoff value
param_OOB = @(x,y) cellfun(@(z)  z.(x),result_data(:,2)) <= y;

EI_valid = param_OOB('EtoI',EI_max);
IE_valid = param_OOB('ItoE',IE_max);
fprintf('\n---E-to-I cutoff = %.2f\n',EI_max)
fprintf('---%i / %i jobs above cutoff\n',sum(~EI_valid),num_files)
fprintf('\n---I-to-E cutoff = %.2f\n',IE_max)
fprintf('---%i / %i jobs above cutoff\n',sum(~IE_valid),num_files)

Tinvalid = Tinvalid | ~EI_valid | ~IE_valid;

fprintf('\n---%i / %i jobs invalid\n',sum(Tinvalid),num_files)
fprintf('\nnum valid jobs = %i\n',sum(~Tinvalid))
fprintf('\nnum jobs to move = %i\n',sum(Tinvalid))


%now move invalid results
f2move = output_fns(Tinvalid);
fprintf('\nmoving files... \n')
%now loop through files & move if invalid
for idx = 1:sum(Tinvalid)
    
    FN = f2move{idx};
    curr_file = fullfile(resdir,FN);
    targ = fullfile(inval_dir,FN);
    movefile(curr_file,targ)
end





