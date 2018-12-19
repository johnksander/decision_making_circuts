clear 
clc
format compact 




%my model 
%---setup---------------------
%jID = 1; %get this from slurm equiv: str2num(getenv('SGE_TASK_ID'));
jID = str2num(getenv('SLURM_ARRAY_TASK_ID'));
t = 5000; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','hpc',...
    'sim_name','parsweep_fastD_baseline','jobID',jID,'tmax',t,...
    'sample_Estay_offset',0,'stim_pulse',[t,0],'cut_leave_state',t);
%---stimulus-----------------
options.stim_targs = 'baseline'; % Eswitch | Estay | baseline
Rstim = 0; %stimulus (in Hz) to target cells 
options.trial_stimuli = [Rstim,Rstim];
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
