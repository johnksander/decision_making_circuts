function Terr = stim_search_wrapper(Tobj,Rstim,options)
%wrapper for fminsearch, used for finding stimulus intensity needed to
%reach a certain average state duration. 
%---takes a target value for mean stimulus duration, stim rate, options
%---returns the error for target duration  

job_runtime = 4; %how long each job will run for, in hours 
Njobs = 1250; %how many jobs to spawn per batch
work = 'run'; %run | debug   just makes commenting/uncommenting stuff less annoying 

options.trial_stimuli = [Rstim,Rstim];%rate for stimulus input spikes

%clear out previous results
switch work
    case 'run'
        update_logfile('cleaning up files...',options.output_log)
        system(sprintf('rm %s/*',options.batchdir));
end


%create a driver file with the current specifications
print_driver_file(options)

partitions = {'paul','guest'};
Nparts = numel(partitions);
Jfile = cell(Nparts,1); %batchfiles 
arrayID = NaN(Nparts,1); %for the job IDs

for idx = 1:Nparts
    %create a SLURM jobfile...
    print_SLURM_file(options,job_runtime,partitions{idx})
    Jfile{idx} = fullfile(options.batchdir,sprintf('batch_job_%s.sh',partitions{idx}));

    %make that file executable
    system(sprintf('chmod +x %s',Jfile{idx}));
end

    
%submit the job
update_logfile('--------------------------',options.output_log)
message = sprintf('---Submitting %i jobs with stimuli = %.1f Hz',Njobs * Nparts,Rstim);
update_logfile(message,options.output_log)
update_logfile(sprintf('---Job durations = %.2f hours...\n',job_runtime),options.output_log)

for idx = 1:Nparts
    switch work
        case 'run'
            [~,thisID] = system(sprintf('sbatch --array=1-%i %s',Njobs,Jfile{idx}));
        case 'debug'
            thisID = 'Submitted batch job 516198 '; warning('work set to debug!');
    end
    
    update_logfile(sprintf('%s job submitted with output: %s',partitions{idx},thisID),options.output_log)
    
    %get the array ID
    thisID = strsplit(thisID);
    thisID = thisID(~cellfun(@isempty,thisID));
    thisID = str2num(thisID{end});
    arrayID(idx) = thisID;
    update_logfile(sprintf('---task ID retained as: %i',arrayID(idx)),options.output_log)
end


switch work
    case 'run'
        %wait for specified amount of time
        wait_time = job_runtime * 3600; % in seconds
        wait_time = wait_time + (2.5 * 60);%add another 2.5 min for good measure
        pause(wait_time)
end

%cancel task ID
update_logfile('stoping remaining jobs...',options.output_log)

for idx = 1:Nparts
    switch work
        case 'run'
            [~,cmd_out] = system(sprintf('scancel %i',arrayID(idx)));
        case 'debug'
            cmd_out = 'job stopped';
    end
    
    update_logfile(sprintf('%s jobs stoped with output: %s',partitions{idx},cmd_out),options.output_log)
end



%error('It''s all good, stopping here')

%submit the job
update_logfile('\n',options.output_log)
update_logfile('::::::: RESULTS :::::::',options.output_log)

%collect jobs results and report 
FNs = dir(fullfile(options.batchdir,'*mat'));
FNs = {FNs.name};
Nfiles = numel(FNs);

message = sprintf('---files found = %i',Nfiles);
update_logfile(message,options.output_log)

state_durations = cell(Nfiles,1);
for idx = 1:Nfiles
    
    data = load(fullfile(options.batchdir,FNs{idx}));
    data = data.sim_results;
    data = data{1};
    valid_states = startsWith(data(:,end),'stim');
    if sum(valid_states) > 0
        data = data(valid_states,2);
        data = cat(1,data{:});
        state_durations{idx} = data;
    end
end

state_durations = cat(1,state_durations{:});
state_durations = state_durations * options.timestep;
Nstates = numel(state_durations);
mu_dur = mean(state_durations);

if Nstates < 1
    update_logfile('',options.output_log)
    update_logfile('WARN!!! NO STATE DURATIONS FOUND',options.output_log)
    update_logfile('---mean stay duration set to = 500s',options.output_log)
    update_logfile('',options.output_log)
    mu_dur = 500;
end



Terr = (Tobj - mu_dur).^2;

message = sprintf('---N states = %i',Nstates);
update_logfile(message,options.output_log)
%print results nicely 
update_logfile('   |   stimulus   |   mean stay   |   error',options.output_log)
message = sprintf('   |    %.1f Hz   |    %.2fs     |   %.2fs (for %.2fs target)',Rstim,mu_dur,Terr,Tobj);
update_logfile(message,options.output_log)
update_logfile('',options.output_log)

% message = sprintf('---trial stimuli = %.1f Hz',Rstim);
% update_logfile(message,options.output_log)
% message = sprintf('---mean stay duration = %.2fs',mu_dur);
% update_logfile(message,options.output_log)
% message = sprintf('---error = %.2fs (for %.2fs target)\n',Terr,Tobj);
% update_logfile(message,options.output_log)




%get the jobID from 

%squeue -u jksander -n fastD_PS -o %F


%squeue -u jksander -n fastD_PS

