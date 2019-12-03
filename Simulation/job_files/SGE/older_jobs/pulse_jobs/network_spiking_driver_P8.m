clear 
clc
format compact 



%my model 
%---setup---------------------
jID = str2num(getenv('SGE_TASK_ID'));
options = set_options('modeltype','PS_stim','comp_location','hpc',...
    'sim_name','sim_v2_P8_1','stim_pulse',[8,1],...
    'jobID',jID,'tmax',755);
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
delete(options.output_log) %no need for these right now