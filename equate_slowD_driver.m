clear
clc
format compact


Tobj = 7.5;
toler = .25; %give it 250 ms tolerance
R0_stim = 25; %start search at 25 hz
search_opt = optimset('TolX',toler);

%my model
%---setup---------------------
num_nets = 10; %number of network pairs
for idx = 1:num_nets %use this to index the different network types
    t = 10e3; %trial simulation time (s)
    options = set_options('modeltype','equate_stim','comp_location','woodstock',...
        'sim_name','equate_slowD_stims','jobID',idx,'tmax',t,...
        'percent_Dslow',.5,'netpair_file','slowD',...
        'stim_pulse',[t,0],'cut_leave_state',t,'sample_Estay_offset',0);
    %---run-----------------------
    fminsearch(@(x) stim_search_wrapper(Tobj,x,options)  ,R0_stim,search_opt)
    
    %Terr = stim_search_wrapper(Tobj,Rstim,options)
    
end

%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,'spikeout_model.m')
end