function [out,info] = find_stay_durations(durations,options,varargin)
%this is meant for analyzing the new (as of 08/27/2018) durations
%datastructure. Takes durations, finds total stay-state durations with
%sampling & "undecided" etc behavior.
%if optional argument 'verify' is given, returns stay-states with valid
%transitions to leave-states. This is meant for supplementing the spikerate
%analysis, since the spike-rate recording conditions were "liberal" (this
%will be fixed for the next job) 

%rules:
%---leave-states must be preceeded by stay-states
%---leave-states must be preceeded by stay-states within the same sample 
%---spontaneous leave-states at the start of a sample are disregarded

if numel(varargin) == 1 && strcmpi(varargin{:},'verify')
    verify_mode = 'on';
else
    verify_mode = 'off';
end


stimulus_duration = options.stim_pulse(1); 
timecourse = size(durations);
timecourse(2) = timecourse(2) + 1; 
timecourse = cell(timecourse);
timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
timecourse(:,3) = cellfun(@(x) mod(x,sum(options.stim_pulse)),timecourse(:,3),'UniformOutput',false);
%current sample's onset, rounding is needed for subsequent operations 
timecourse(:,4) = cellfun(@(x,y) round(x-y,2),timecourse(:,1),timecourse(:,3),'UniformOutput',false);
samp_onsets = unique(cat(1,timecourse{:,4})); %like unique won't work properly here without rounding 
timecourse(:,4) = cellfun(@(x) find(x==samp_onsets),timecourse(:,4),'UniformOutput',false); %would also break without rounding
timecourse(:,end) = durations(:,end);
%timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
%timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});
leave_states = find(strcmp(timecourse(:,end),'leave'));
invalid_switches = zeros(size(leave_states));
for idx = 1:numel(leave_states)
    %leave_durr = timecourse{leave_states(idx),2}; %How long leave-state lasted
    %leave_start = timecourse{leave_states(idx),1} - leave_durr; %when that leave state started
    leave_sample = timecourse{leave_states(idx),4};
    curr_rec = timecourse(1:leave_states(idx),:); %duration record up to that point
    %check for no prior stay-state (2 sequential leaves after 1st artificial stay)
    %or for leave states prior to last stay (sequential leaves during normal sim)
    prev_stay = find(startsWith(curr_rec(:,end),'stim'), 1, 'last');
    prev_leave = [0;leave_states]; %just so this wont break with idx == 1
    prev_leave = max(prev_leave(prev_leave < leave_states(idx)));
    if isempty(prev_stay) || prev_leave > prev_stay
        invalid_switches(idx) = 1;
    else
        prev_stay = curr_rec(prev_stay,:); %stay-state prior to leave-transition
        prev_stay_sample = prev_stay{4};
        if leave_sample ==  prev_stay_sample 
            %a stay-state existed in this sample before transitioning to leave
            %i.e. legitimate switch in the middle of a sample
        else
            %it's bogus or something
            invalid_switches(idx) = 1;
        end
    end
end

invalid_switches = leave_states(invalid_switches == 1);
timecourse(invalid_switches,:) = []; %remove invalid leaves
undecided_states = strcmp(timecourse(:,end),'undecided');
timecourse(undecided_states,:) = []; %remove undecided states 
%ok... now for the real results 
leave_states = find(strcmp(timecourse(:,end),'leave'));
out = NaN(numel(leave_states),3);
out = array2table(out,'VariableNames',{'duration','decision_time','samples'});
switch verify_mode
    case 'on'
        %{state.timeidx,state.count,state.sample_clock,record_stim_label});
        info = cell(numel(leave_states),2); %just return the time index & label for verification
end
leave_states = cat(1,0,leave_states); %just for indexing 
for idx = 2:numel(leave_states) %start at 2 b/c you appended that zero for indexing 
    curr_rec = leave_states(idx-1)+1:leave_states(idx);
    curr_rec = timecourse(curr_rec,:); %stay states up to that point 
    total_duration = curr_rec(1:end-1,2); %durations of all the stay-states (not inc. leave...)
    total_duration = sum(cat(1,total_duration{:}));
    out.duration(idx-1) = total_duration;
    %stay_start = curr_rec{end-1,1} - curr_rec{end-1,2}; %when previous stay started 
    %leave_start = curr_rec{end,1} - curr_rec{end,2}; %when this leave state started (aka decision time)
    leave_dur = curr_rec{end,2};
    dec_time = curr_rec{end,3}; %use the sample clock for this, more accurate if there's  
    %multiple "samples" within the same pulse i.e. stay-undecided-stay (undecided is usally extremely brief, probably
    %just momentarily below our threshold) 
    out.decision_time(idx-1) = dec_time - leave_dur; 
    total_samples = cat(1,curr_rec{:,4}); %use the unique sample numbers
    total_samples = numel(unique(total_samples));
    out.samples(idx-1) = total_samples;
    switch verify_mode
        case 'on'
            info(idx-1,:) = curr_rec(end-1,[1,5]);%previous state's time index & label for verification
    end
end











   