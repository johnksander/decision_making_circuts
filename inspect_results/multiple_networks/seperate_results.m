clear;clc
format compact
close all

%Do you want to seperate simulation results according to some criteria?
%Then this is the place for you. 



%this must be parfored... too many files

basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

num_workers = 24;
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf) %,'AttachedFiles',{which('find_stay_durations')})



sim_name = 'nets_D2t_pref';
seperate_into = 'stimB-%i'; %B is what percent of A

resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,'*.mat')); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);

special_progress_tracker = fullfile(basedir,'inspect_results','multiple_networks','SPT.txt');
if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start

num_files = numel(output_fns);
fprintf('-----starting simulation %s (%i files)\n',sim_name,num_files)

parfor idx = 1:num_files
    
    FN = output_fns{idx};
    curr_file = load(FN);
    
    %----seperation criteria
    stims = curr_file.options.trial_stimuli;
    B = stims(2) ./stims(1);
    %where it goes
    if round(B*100) == 100
        outdir = sprintf('%s_%s',resdir,sprintf(seperate_into,round(B*100)));
        if ~isdir(outdir),mkdir(outdir);end
        copyfile(FN,outdir)
    end
    
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .1)) == 0 %at ten percent
        progress = (progress / num_files) * 100;
        fprintf('%s ---- %.1f percent complete\n',datestr(now,31),progress);
    end
    
end

delete(special_progress_tracker)

delete(gcp('nocreate'))





