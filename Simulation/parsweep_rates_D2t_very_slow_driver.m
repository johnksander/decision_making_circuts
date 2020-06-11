clear 
clc
format compact 

%get the rates for these networks

sim_name = 'parsweep_D2t-slower_spikerates'; %annoying but needed...

%my model 
%---setup---------------------
jID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
t = 25; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','hpc',...
    'sim_name',sim_name,'jobID',jID,'tmax',t,...
    'ratelim_check','on','cut_leave_state',t,'noswitch_timeout',t+1);

%files should already be here
FN = dir(fullfile(options.save_dir,'*mat'));
FN = {FN.name};
FN = FN{jID};
FN = fullfile(options.save_dir,FN);

conn_params = load(FN);
conn_params = conn_params.options;
options.EtoI = conn_params.EtoI;
options.ItoE = conn_params.ItoE;
options.jobID = conn_params.jobID; %VERY IMPORTANT
options.sim_name = sprintf('PS_%s_%i',sim_name,options.jobID);

%---run-----------------------
run_job = true;
while run_job

    modelfile = special_ratejob_model(options);
    results = load(FN);
    if isfield(results,'sim_results')
        results = results.sim_results;
        results = results(4); %should be where ratelim is 
        run_job = isempty(results); %stop running if you've got rate data 
    end
end


%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
delete(options.output_log) %no need for these right now
% logdir = fullfile(options.save_dir,'logs'); %put them seperately
% if ~isdir(logdir),mkdir(logdir);end
% movefile(options.output_log,logdir)