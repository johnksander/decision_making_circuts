clear
clc
format compact


Tobj = 7.5;
toler = .25; %give it 250 ms tolerance
R0_stim = 25; %start search at 25 hz
search_opt = optimset('TolX',toler);t = 1; %trial simulation time (s)
num_nets = 10; %number of network pairs


%---setup---------------------
for idx = 1:num_nets %use this to index the different network types
    
    %:::start:::
    t = 50; %trial simulation time (s) 
    options = set_options('modeltype','equate_stim','comp_location','hpc',...
        'sim_name','equate_slowD_stims','jobID',idx,'tmax',t,...
        'percent_Dslow',.5,'netpair_file','slowD',...
        'stim_pulse',[t,0],'sample_Estay_offset',0);
    %:::end:::
    
    %check if network has been optimized yet
    FN = fullfile(options.save_dir,sprintf('%s.mat',options.sim_name));
    if exist(FN,'file') == 0 %this network has not been optimized
        options.master_driver = which(mfilename);
        
        %---run-----------------------
        fminsearch(@(x) stim_search_wrapper(Tobj,x,options)  ,R0_stim,search_opt)
        
        %Terr = stim_search_wrapper(Tobj,Rstim,options)
        
        
        %---cleanup-------------------
        if isempty(dir(fullfile(options.save_dir,'code4*zip')))
            driverfile = mfilename;
            backup_jobcode(options,driverfile,'spikeout_model.m')
        end        
    end
end

%note---
% the start & end comments are used by stim_search_wrapper to
% identify what parameters must be specified in the driver file. 
% there can be no comments in that section!!
