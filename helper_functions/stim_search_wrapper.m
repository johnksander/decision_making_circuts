function Terr = stim_search_wrapper(Tobj,Rstim,options)
%wrapper for fminsearch, used for finding stimulus intensity needed to
%reach a certain average state duration. 
%---takes a target value for mean stimulus duration, stim rate, options
%---returns the error for target duration  


options.trial_stimuli = [Rstim,Rstim];%rate for stimulus input spikes

%---run-----------------------
spikeout_model(options);

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
