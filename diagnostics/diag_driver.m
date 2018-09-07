clear
clc
format compact
hold off;close all
%investigating model behavior

addpath('../')
%jobID = str2num(getenv('SGE_TASK_ID'));

%---setup---------------------
t = 50;
options = set_options('modeltype','PS_stim',...
    'sim_name','diagnostics_Sgtest_transients',...
    'jobID',5,'tmax',t,'stim_pulse',[t,0],'stim_schedule','flexible',...
    'comp_location','woodstock','cut_leave_state',5e-3);
stim_mod = .85:.03:1.15; 
%just randomly pick one from the spread, do enough it'll even out 
stim_mod = randsample(stim_mod,1);
%adjust stim B
options.trial_stimuli(2) = options.trial_stimuli(2) * stim_mod; 
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model(options);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(options,driverfile,modelfile)
delete(options.output_log) %no need for these right now



