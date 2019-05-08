function [out,info] = find_stay_durations(durations,options,varargin)
%this is meant for analyzing the new (as of 08/27/2018) durations
%datastructure. Takes durations, finds total stay-state durations with
%sampling & "undecided" etc behavior.
%if optional argument 'verify' is given, returns stay-states with valid
%transitions to leave-states. This was meant for supplementing the spikerate
%analysis, since the spike-rate recording conditions was "liberal". 
%Durations input format: {timestamp, state count, sample clock, stim label}

%rules:
%---leave-states must be preceeded by stay-states
%---leave-states must be preceeded by stay-states within the same sample 
%---spontaneous leave-states at the start of a sample are disregarded

if numel(varargin) == 1 && strcmpi(varargin{:},'verify')
    verify_mode = 'on';
else
    verify_mode = 'off';
end


if all(isnan(options.stim_pulse)) %there's a legit pulsing stim 
    sample_duration = options.tmax;
else %constant stim 
    sample_duration = sum(options.stim_pulse);
end

timecourse = size(durations);
timecourse(2) = timecourse(2) + 1; 
timecourse = cell(timecourse);
timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
timecourse(:,3) = cellfun(@(x) mod(x,sample_duration),timecourse(:,3),'UniformOutput',false);
%must use uniquetol() and ismembertol() for roundoff errors 
timecourse(:,4) = cellfun(@(x,y) x-y,timecourse(:,1),timecourse(:,3),'UniformOutput',false);
samp_onsets = uniquetol(cat(1,timecourse{:,4})); 
timecourse(:,4) = cellfun(@(x) find(ismembertol(samp_onsets,x)),timecourse(:,4),'UniformOutput',false); 
timecourse(:,end) = durations(:,end);
timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});
leave_states = find(strcmp(timecourse.state,'leave'));
invalid_switches = zeros(size(leave_states));
for idx = 1:numel(leave_states)
    %leave_durr = timecourse{leave_states(idx),2}; %How long leave-state lasted
    %leave_start = timecourse{leave_states(idx),1} - leave_durr; %when that leave state started
    leave_sample =  timecourse.sample_number(leave_states(idx));
    curr_rec = timecourse(1:leave_states(idx),:); %duration record up to that point
    %check for no prior stay-state (2 sequential leaves after 1st artificial stay)
    %or for leave states prior to last stay (sequential leaves during normal sim)
    prev_stay = find(startsWith(curr_rec.state,'stim'), 1, 'last');
    prev_leave = [0;leave_states]; %just so this wont break with idx == 1
    prev_leave = max(prev_leave(prev_leave < leave_states(idx)));
    if isempty(prev_stay) || prev_leave > prev_stay
        invalid_switches(idx) = 1;
    else
        prev_stay = curr_rec(prev_stay,:); %stay-state prior to leave-transition
        prev_stay_sample = prev_stay.sample_number;
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
undecided_states = strcmp(timecourse.state,'undecided');
timecourse(undecided_states,:) = []; %remove undecided states 
%ok... now for the real results 
leave_states = find(strcmp(timecourse.state,'leave'));
out = NaN(numel(leave_states),3);
out = array2table(out,'VariableNames',{'duration','decision_time','samples'});
switch verify_mode
    case 'on'
        %{state.timeidx,state.count,state.sample_clock,record_stim_label});
        info = cell2table(cell(numel(leave_states),3),'VariableNames',{'event_time','state','spiking_data'});
        %return the time index, label, whether there's good spiking data
end

leave_states = cat(1,0,leave_states); %just for indexing 
for idx = 2:numel(leave_states) %start at 2 b/c you appended that zero for indexing 
    curr_rec = leave_states(idx-1)+1:leave_states(idx);
    curr_rec = timecourse(curr_rec,:); %stay states up to that point 
    total_duration = curr_rec.duration(1:end-1); %durations of all the stay-states (not inc. leave...)
    out.duration(idx-1) = sum(total_duration);
    %stay_start = curr_rec{end-1,1} - curr_rec{end-1,2}; %when previous stay started 
    %leave_start = curr_rec{end,1} - curr_rec{end,2}; %when this leave state started 
    curr_leave = curr_rec(end,:);
    dec_time = curr_leave.sample_time; %use the sample clock for this, more accurate if there's  
    %multiple "samples" within the same pulse i.e. stay-undecided-stay (undecided is usally extremely brief, probably
    %just momentarily below our threshold) 
    out.decision_time(idx-1) = dec_time - curr_leave.duration; 
    total_samples = numel(unique(curr_rec.sample_number));%use the unique sample numbers
    out.samples(idx-1) = total_samples;
    switch verify_mode
        case 'on'
            info.event_time{idx-1} = curr_rec.event_time(end-1);
            info.state{idx-1} = curr_rec.state{end-1};
            
            %---check if spiking data is within valid recording window
            prev_stay = curr_rec(end-1,:); %stay-state prior to leave-transition
            leave_start = curr_leave.event_time - curr_leave.duration;
            %stay-state followed by a leave-state within Yms later, and lasted at least Y ms
            recwin_check = leave_start - prev_stay.event_time;
            recwin_check = recwin_check <= options.record_postswitch ...
                &&  (recwin_check + curr_leave.duration) >= options.record_postswitch;
            %if the state-state lasted at least Xms, and meets above criteria
            recwin_check = prev_stay.duration > options.record_preswitch && recwin_check;
            info.spiking_data{idx-1} = recwin_check;
    end
end











   