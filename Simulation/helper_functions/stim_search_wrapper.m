function Terr = stim_search_wrapper(Tobj,Rstim,options)
%wrapper for fminsearch, used for finding stimulus intensity needed to
%reach a certain average state duration. 
%---takes a target value for mean stimulus duration, stim rate, options
%---returns the error for target duration  

job_runtime = .75; %how long each job will run for, in hours 
Njobs = 1e3; %how many jobs to spawn per batch
work = 'run'; %run | debug   just makes commenting/uncommenting stuff less annoying 


options.trial_stimuli = [Rstim,Rstim];%rate for stimulus input spikes

%create a driver file with the current specifications
print_driver_file(options)

%create a SLURM jobfile...
print_SLURM_file(options,job_runtime)
Jfile = fullfile(options.batchdir,'batch_job.sh');


%make that file executable
system(sprintf('chmod +x %s',Jfile));

%submit the job
update_logfile('--------------------------',options.output_log)
message = sprintf('---Submitting %i jobs with stimuli = %.1f Hz',Njobs,Rstim);
update_logfile(message,options.output_log)
update_logfile(sprintf('---Job durations = %.1f hours...\n',job_runtime),options.output_log)

switch work
    case 'run'
        [~,arrayID] = system(sprintf('sbatch --array=1-%i %s',Njobs,Jfile));
    case 'debug'
        arrayID = 'Submitted batch job 516198';
end

update_logfile(sprintf('job submitted with output: %s',arrayID),options.output_log)

%get the array ID 
arrayID = strsplit(arrayID);
arrayID = arrayID(~cellfun(@isempty,arrayID));
arrayID = str2num(arrayID{end}); 
update_logfile(sprintf('task ID number retained as: %i',arrayID),options.output_log)


%wait for specified amount of time 
wait_time = job_runtime * 3600; % in seconds 
wait_time = wait_time + (5 * 60);%add another 5 min for good measure
pause(wait_time)


%cancel task ID
update_logfile('stoping remaining jobs...',options.output_log)

switch work
    case 'run'
        [~,cmd_out] = system(sprintf('scancel %i',arrayID));
    case 'debug'
        cmd_out = 'job stopped';
end

update_logfile(sprintf('jobs stoped with with output: %s',cmd_out),options.output_log)


error('It''s all good, stopping here')
%collect jobs results and report 


%get the jobID from 

%squeue -u jksander -n fastD_PS -o %F


%squeue -u jksander -n fastD_PS


%notes: 
%set error to high if no datapoints, etc 
%paul has a 3 day timelimit

%---run-----------------------
%spikeout_model(options);

state_durations = fullfile(options.save_dir,options.sim_name);
state_durations = load(state_durations);
state_durations = state_durations.sim_results;
state_durations = state_durations{1};
valid_states = startsWith(state_durations(:,end),'stim');
state_durations = state_durations(valid_states,2);
state_durations = cat(1,state_durations{:});
state_durations = state_durations * options.timestep;
mu_dur = mean(state_durations);
Terr = Tobj - mu_dur;

message = sprintf('---trial stimuli = %.1f Hz',Rstim);
update_logfile(message,options.output_log)
message = sprintf('---mean stay duration = %.2fs',mu_dur);
update_logfile(message,options.output_log)
message = sprintf('---error = %.2fs (for %.2fs target)\n',Terr,Tobj);
update_logfile(message,options.output_log)
