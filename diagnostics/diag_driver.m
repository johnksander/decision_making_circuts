clear
clc
format compact

%investigating model behavior

addpath('../')
%jobID = str2num(getenv('SGE_TASK_ID'));

%---setup---------------------
options = set_options('modeltype','PS_stim',...
    'sim_name','diagnostics_Sgtest',...
    'jobID',1,'tmax',20,'stim_pulse',[1,1],...
    'force_back2stay','true',...
    'comp_location','woodstock');


%NOTE: I AM NOT removing the first artificial stay state here. 
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



