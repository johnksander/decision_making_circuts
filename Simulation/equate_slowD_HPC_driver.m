clear
clc
format compact

%NOTE: using longer dt for this search 
Tobj = 7.5; %target mean duration 
R0_stim = 25; %start search at 25 hz

%stopping criteria (both must be met)
fun_tol = .25 .^2; % 250 ms tolerance for changes in objective function
X_tol = 2; % 2 Hz tolerance for change in stimulus (per step)
search_opt = optimset('TolFun',fun_tol,'TolX',X_tol);

num_nets = 10; %number of network pairs

%---setup---------------------
for idx = 1:num_nets %use this to index the different network types
    
    %:::start:::
    t = 200; %trial simulation time (s)
    options = set_options('modeltype','equate_stim','comp_location','hpc',...
        'timestep',.25e-3,...
        'sim_name','equate_slowD_stims','jobID',idx,'tmax',t,...
        'percent_Dslow',.5,'netpair_file','slowD',...
        'stim_pulse',[t,0],'sample_Estay_offset',0);
    %:::end:::
    
    %check if network has been optimized yet
    FN = fullfile(options.save_dir,sprintf('%s.mat',options.sim_name));
    if exist(FN,'file') == 0 %this network has not been optimized
        
        options.master_driver = which(mfilename);
        solution = false;
        
        while ~solution
            %---run-----------------------
            
            [Req,~,exitflag] = ...
                fminsearch(@(x) stim_search_wrapper(Tobj,x,options)  ,R0_stim,search_opt);
            
            solution = exitflag == 1; %ensure minima actually found
        end
        
        %solution found, collect results and delete batchfiles
        %submit the job
        update_logfile('********************************',options.output_log)
        update_logfile('********Soulution found*********',options.output_log)
        update_logfile('',options.output_log)
        
        %collect jobs results and report
        bfiles = dir(fullfile(options.batchdir,'*mat'));
        bfiles = {bfiles.name};
        Nfiles = numel(bfiles);
        state_durations = cell(Nfiles,1);
        for fileidx = 1:Nfiles
            data = load(fullfile(options.batchdir,bfiles{fileidx}));
            data = data.sim_results;
            data = data{1};
            valid_states = startsWith(data(:,end),'stim');
            if sum(valid_states) > 0
                data = data(valid_states,2);
                data = cat(1,data{:});
                state_durations{fileidx} = data;
            end
        end
        state_durations = cat(1,state_durations{:});
        state_durations = state_durations * options.timestep;
        %save durations, equated stim value, options
        save(FN,'state_durations','Req','options')
        %remove batch files
        system(sprintf('rm -r %s',options.batchdir));
        update_logfile('file cleanup complete',options.output_log)
    end
end

%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,'spikeout_model.m')
end

%note---
% the start & end comments are used by stim_search_wrapper to
% identify what parameters must be specified in the driver file.
% there can be no comments in that section!!

% %use this for like... debugging n crap
% [Req,~,exitflag] = ...
%     fminsearch(@(x) ((x-100).^2) * all([x < 500,x > -500]) ,R0_stim);
