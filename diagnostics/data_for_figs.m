clear
clc
format compact
hold off;close all
%show network characteristics


addpath('../')
Sname = 'example_behavior_equated';

jobs = [3:13,4:14]; %do a few runs for slow net #2

for idx = 1:numel(jobs)
    
    jobID = jobs(idx);
    
    %---setup---------------------
    tmax = 60;
    options = set_options('modeltype','diagnostics','comp_location','woodstock',...
        'sim_name',Sname,'jobID',jobID,'tmax',tmax,'netpair_file','D2t-slower',...
        'noswitch_timeout',tmax,'cut_leave_state',tmax+1);
    
    do_net = mod(options.jobID,10);
    do_net(do_net == 0) = 10;
    options = get_network_params(do_net,options);
    options.EtoE = .0405; %fixed
    %options.trial_stimuli{1} = [0,0]; %no stimulus
    
    %---run-----------------------
    exit_status = false;
    while ~exit_status
        [modelfile,exit_status] = diag_model_lite(options);
    end
    %---cleanup-------------------
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
    delete(options.output_log) %no need for these right now
    
    fprintf('skipping inspect() call!\n')
    %setenv('JID',num2str(options.jobID));
    %setenv('SIM_NAME',Sname);
    %inspect
    
    hold off;close all

    fprintf('job finished (JID = %i)\n',options.jobID)
end
