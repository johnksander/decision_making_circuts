clear 
clc
format compact 




%my model 
%---setup---------------------
jID = str2num([getenv('SLURM_JOBID'), getenv('SLURM_ARRAY_TASK_ID')]);
t = 1500; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','hpc',...
    'sim_name','parsweep_D2t_very_slow_baseline','jobID',jID,'tmax',t,...
    'ratelim_check','off','cut_leave_state',t,'noswitch_timeout',t);


%---grid search
Ngrid = 100;
ItoE = linspace(0.1,12.5,Ngrid);
EtoI = linspace(0,.75,Ngrid);
[ItoE,EtoI] = meshgrid(ItoE,EtoI);
ItoE = ItoE(:); EtoI = EtoI(:);
%looking only for slow networks
valid_range = ItoE >= 7.5 & EtoI >= .225;
ItoE = ItoE(valid_range);
EtoI = EtoI(valid_range);
options.grid_index = str2num(getenv('SLURM_ARRAY_TASK_ID'));%HPCC only lets indicies up to 10k!!
options.ItoE = ItoE(options.grid_index);
options.EtoI = EtoI(options.grid_index);

%---run-----------------------
modelfile = spikeout_model_grid(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)