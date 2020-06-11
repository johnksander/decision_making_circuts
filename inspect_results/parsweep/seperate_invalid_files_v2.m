clear
clc
format compact
hold off;close all


%This code is intended to move result files that do not match some criteria
%(e.g. state duration time, parameter values, etc) into a subdirectory.
%This preserves all results from a parameter sweep, while simplifying
%analyses after determining those criteria. 

sim_name = 'parsweep_D2t-slower_spikerates'; %MUST CHANGE CRITERIA TOO

%----CRITERIA----
%must be in the summary results file (summary file from parsweep_D2t_very_slow_baseline)
jobs = '/Users/ksander/Desktop/work/ACClab/rotation/project/Results/summary_file.mat';
jobs = load(jobs);
jobs = jobs.result_data;
%jobs = jobs(:,2);


%specify simulation
%---sim setup-----------------
basedir = '~/Desktop/work/ACClab/rotation/project'; %'/home/acclab/Desktop/ksander/rotation/project';
resdir = fullfile(basedir,'Results',sim_name);
finished_dir = fullfile(resdir,'already_finished');
if ~isdir(finished_dir),mkdir(finished_dir);end
%addpath(fullfile(basedir,'helper_functions'))

t = 25; %trial simulation time (s) 
%options = set_options('modeltype','PS','comp_location','hpc',...
%   'sim_name','parsweep_D2t-slower_spikerates','jobID',jID,'tmax',t,...
%    'ratelim_check','on','cut_leave_state',t,'noswitch_timeout',t);




num_files = numel(jobs(:,1));
for idx = 1:num_files
    curr_file = jobs(idx,:);
    has_rates = ~cellfun(@isempty,curr_file(3));
    FN = sprintf('PS_parsweep_D2t-slower_spikerates_%i.mat',idx);
    makefile = true;
    if exist(fullfile(resdir,FN)) > 0
        Fdata = load(fullfile(resdir,FN));
        if isfield(Fdata,'sim_results')
            Fdata = Fdata.sim_results;
            Fdata = Fdata(4); %should be where ratelim ris
            if ~isempty(Fdata)
                 movefile(fullfile(resdir,FN),finished_dir)
                 makefile = false;
                 disp('found one')
            end
        end
    end
    
    sim_results = cell(1,4);
    options = curr_file{2};
    options.jobID = idx;
    options.tmax = t;
    options.cut_leave_state = t;
    options.noswitch_timeout = t;
    options.sim_name = 'parsweep_D2t-slower_spikerates';
    options.ratelim.check = 'on';
    if has_rates
        sim_results{4} = curr_file{3};
        save(fullfile(finished_dir,FN),'options','sim_results')
    elseif makefile
        save(fullfile(resdir,FN),'options')
        
    elseif has_rates && ~makefile
        error('what do I do here')
    end
end






