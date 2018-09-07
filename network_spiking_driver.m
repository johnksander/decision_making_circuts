clear 
clc
format compact 



%my model 
%---setup---------------------
t = 1e3; %1k time 
jID = str2num(getenv('SGE_TASK_ID'));
options = set_options('modeltype','PS_stim','comp_location','hpc',...
    'sim_name','dep_pref','record_spiking','off',...
    'jobID',jID,'tmax',t,'stim_pulse',[t,0],...
    'stim_schedule','flexible','cut_leave_state',5e-3);
%adjust stimulus B strength
stim_mod = .85:.03:1.15; %just randomly sample mod weight, do enough it'll even out 
stim_mod = randsample(stim_mod,1);
options.trial_stimuli(2) = options.trial_stimuli(2) * stim_mod; %adjust stim B
%---run-----------------------
modelfile = spikeout_model(options);
%---cleanup-------------------
if options.jobID <= 10 %only do this for one set...
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
delete(options.output_log) %no need for these right now