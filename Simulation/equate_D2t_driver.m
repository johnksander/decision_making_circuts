clear
clc
format compact

%NOTE: using longer dt for this search
Tobj = 7.5; %target mean duration

%target time tolerance (must pass check to save)
targ_tol = .11 ^2; %110 ms tolerance

%stopping criteria (both must be met)
%fun_tol = .25 ^2; % 250 ms tolerance for changes in objective function
%X_tol = .25; % .2 Hz tolerance for change in stimulus (per step)
%search_opt = optimset('TolFun',fun_tol,'TolX',X_tol); fminbnd()

num_nets = 10; %number of network pairs

%---setup---------------------
for idx = 1:num_nets %use this to index the different network types
    
    %:::start:::
    t = 350; %trial simulation time (s)
    options = set_options('modeltype','equate_stim','comp_location','hpc',...
        'sim_name','equate_D2t_stims-try2','timestep',.25e-3,...
        'netpair_file','D2t','jobID',idx,'tmax',t);
    %:::end:::
    
    %check if network has been optimized yet
    FN = fullfile(options.save_dir,sprintf('%s.mat',options.sim_name));
    if exist(FN,'file') == 0 %this network has not been optimized
        
        options.master_driver = which(mfilename);
        stop_search = false;
        solution = false;
        
        %set bounds based on previous search results
        Rprev = unique(options.trial_stimuli);
        Rmax = Rprev + .5*Rprev;
        Rmin = Rprev - .5*Rprev;
        
        
        srgFN = fullfile(options.save_dir,sprintf('%s_checkpoint.mat',options.sim_name));
        search_opt = optimoptions('surrogateopt','CheckpointFile',srgFN,...
            'Display','iter','PlotFcn',[],'InitialPoints',Rprev);
        if exist(search_opt.CheckpointFile) > 0
            fprintf('\n RESUMING surrogate search from file:\n::: %s\n',search_opt.CheckpointFile)
            start_new = false;
        else
            start_new = true;
        end
        
        
        while ~stop_search
            %---run-----------------------
            
            %either do fminsearch()
            %[Req,~,exitflag] = ...
            %    fminsearch(@(x) stim_search_wrapper(Tobj,x,options)  ,R0_stim,search_opt);
            
            %%or fminbnd()
            %[Req,Terr,exitflag] = ...
            %    fminbnd(@(x) stim_search_wrapper(Tobj,x,options),Rmin,Rmax,search_opt);
            
            
            %surrogateopt()
            if start_new
                
                [Req,Terr,exitflag] = ...
                    surrogateopt(@(x) stim_search_wrapper(Tobj,x,options),Rmin,Rmax,search_opt);
                
            else %resume from checkpoint file
                [Req,Terr,exitflag] = surrogateopt(search_opt.CheckpointFile);
                
            end
            
            stop_search = exitflag == 1; %see if search alg is done
            if Terr < targ_tol
                solution = true; %found a good parameter
            end
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
        if solution %only save the completed datafile if you found a result within toleranace range!
            %save durations, equated stim value, options
            save(FN,'state_durations','Req','options')
        end
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
