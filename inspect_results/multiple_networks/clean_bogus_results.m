clear
clc
format compact

%this must be parfored... too many files

Snames = {'nets_mixstim-NOBSTEST'}; %Snames = {'nets_fastD','nets_slowD'};

basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

num_workers = 24;
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf) %,'AttachedFiles',{which('find_stay_durations')})

for sidx = 1:numel(Snames)
    fix_this(basedir,Snames{sidx})
end

delete(gcp('nocreate'))


function fix_this(basedir,sim_name)
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
num_deleted = false(num_files,1);
parfor idx = 1:num_files
    %if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %get state durations
    try
        state_durations = curr_file.sim_results;
        %--plus these steps now
        state_durations = state_durations{1};
        [state_durations,Sinfo] = find_stay_durations(state_durations,curr_file.options,'verify');
    catch
        
        delete(output_fns{idx})
        num_deleted(idx) = true;
    end
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .05)) == 0 %at half a percent
        progress = (progress / num_files) * 100;
        fprintf('%s ---- %.1f percent complete\n',datestr(now,31),progress);
    end
    
end

fprintf('\n\n-----total deleted: %i out of %i total files\n\n',sum(num_deleted),num_files)
delete(special_progress_tracker)

end
