clear;clc
format compact
close all

%this was intended for jobs that starting running w/ spiking data saved,
%but then stopped saving spiking data (or vice-versa). Finds the results
%with spiking data and sets them aside 



%this must be parfored... too many files

basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

num_workers = 24;
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf) %,'AttachedFiles',{which('find_stay_durations')})


sim_name = 'nets_D2t_pref';




resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_baseline'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);

special_progress_tracker = fullfile(basedir,'inspect_results','multiple_networks','SPT.txt');
if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start

fprintf('-----starting simulation results: %s\n',sim_name)

num_files = numel(output_fns);
num_found = false(num_files,1);
parfor idx = 1:num_files
    %if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %get state durations
    switch curr_file.options.record_spiking
        case 'off'
        case 'on'
             num_found(idx) = true;
    end

    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .05)) == 0 %at half a percent
        progress = (progress / num_files) * 100;
        fprintf('%s ---- %.1f percent complete\n',datestr(now,31),progress);
    end
    
end

fprintf('\n\n-----total found: %i out of %i total files\n\n',sum(num_found),num_files)
delete(special_progress_tracker)

delete(gcp('nocreate'))

targ_dir = [resdir,'_spikedata'];
if ~isdir(targ_dir),mkdir(targ_dir);end
F2copy = output_fns(num_found);

fprintf('copying files...\n')
cellfun(@(x) copyfile(x,targ_dir),F2copy)











