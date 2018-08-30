clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
%jobID = str2num(getenv('SGE_TASK_ID'));

%---setup---------------------
options = set_options('modeltype','PS_stim',...
    'sim_name','diagnostics_Sgtest_transients',...
    'jobID',1,'tmax',50,'stim_pulse',[6,1],'stim_schedule','flexible',...
    'comp_location','woodstock');

%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



