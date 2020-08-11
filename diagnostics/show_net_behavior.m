clear
clc
format compact
hold off;close all
%show network characteristics


addpath('../')
Sname = 'example_behavior';

jobs = 14:10:44; %do a few runs for slow net #2
tmax = 25;

for idx = 1:numel(jobs)
    
    %---setup---------------------
    options = set_options('modeltype','diagnostics','comp_location','woodstock',...
        'sim_name',Sname,'jobID', jobs(idx),'tmax',tmax,'netpair_file','D2t-slower',...
        'noswitch_timeout',tmax+1,'cut_leave_state',tmax);
    
    do_net = mod(options.jobID,10);
    do_net(do_net == 0) = 10;
    options = get_network_params(do_net,options);
    options.EtoE = .0405; %fixed
    options.trial_stimuli{1} = [0,0]; %no stimulus
    
    run_this_job(Sname,options)
    
    fprintf('job finished (JID = %i)\n',options.jobID)
    
end

function run_this_job(Sname,opt)
%---run-----------------------
exit_status = false;
while ~exit_status
    [modelfile,exit_status] = diag_model_lite(opt);
end
%---cleanup-------------------
driverfile = mfilename;
backup_jobcode(opt,driverfile,modelfile)
delete(opt.output_log) %no need for these right now

setenv('JID',num2str(opt.jobID));
setenv('SIM_NAME',Sname);
inspect
hold off;close all
end
