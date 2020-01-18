clear
clc
format compact
hold off;close all


%This code is intended to move result files that do match some criteria
%(e.g. state duration time, parameter values, etc) into a subdirectory.
%This preserves all results from a parameter sweep, while simplifying
%analyses after determining those criteria. 

sim_name = 'parsweep_D2t_baseline'; %MUST CHANGE CRITERIA TOO

%----CRITERIA----
%connection strengths
%EI_max = .75; IE_max = 12.5; %these set connection maximums
Ngrid = 100;
ItoE = linspace(0.1,12.5,Ngrid);
EtoI = linspace(0,.75,Ngrid);
[ItoE,EtoI] = meshgrid(ItoE,EtoI);
ItoE = ItoE(:); EtoI = EtoI(:);
valid_range = ItoE >= 7.5 & EtoI >= .225;%this was "very slow parsweep"
valid_range = ~valid_range; %so get everything we didn't search over in that job 
valid_params = [ItoE(valid_range),EtoI(valid_range)];


num_workers = 24; %for parfor loading results 

%specify simulation
%---sim setup-----------------
basedir = '/home/acclab/Desktop/ksander/rotation/project';
resdir = fullfile(basedir,'Results',sim_name);
merge_dir = fullfile(resdir,'merge_results');
addpath(fullfile(basedir,'helper_functions'))
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = {output_fns.name};
if ~isdir(merge_dir),mkdir(merge_dir);end
%get results
num_files = numel(output_fns);
%parfor stuff
output_log = fullfile(resdir,'output_log.txt');
special_progress_tracker = fullfile(resdir,'SPT.txt');
if exist(special_progress_tracker) > 0
    delete(special_progress_tracker) %fresh start
end
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf)

valid_jobs = false(num_files,1);
parfor idx = 1:num_files
    %if mod(idx,1000) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    FN = fullfile(resdir,output_fns{idx});
    curr_file = load(FN);
    %job params
    params = curr_file.options;
    params = [params.ItoE,params.EtoI];
    if ismember(params,valid_params,'rows')
        %copy to the merge dir
        copyfile(FN,merge_dir)
        valid_jobs(idx) = true;
    end
 
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .1)) == 0 %at 10 percent
        progress = (progress / num_files) * 100;
        message = sprintf('----%.1f percent complete',progress);
        update_logfile(message,output_log)
    end
end
delete(gcp('nocreate'));delete(special_progress_tracker);delete(output_log)

fprintf('\n %i of %i total files copied \n',sum(valid_jobs),num_files)

