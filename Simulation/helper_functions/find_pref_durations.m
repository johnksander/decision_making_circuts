function [valid_data,Tinfo] = find_pref_durations(durations,options,varargin)
%this is meant for analyzing the new (as of 08/27/2018) durations
%datastructure. Takes durations, finds stay-state durations with
%sampling & "undecided" etc behavior.
%Also returns "Tinfo", gives information about trainsition types &
%latencies (i.e. how long between A stops and B begins)


%rules:
%---leave-states must be preceeded by stay-states
%---leave-states must be preceeded by stay-states within the same sample
%---spontaneous leave-states at the start of a sample are disregarded

%add "inter-stimulus-interval" information output here

%for strict rules, stimA-undecided-stimA, only use last stimA state
%for relaxed rules, stimA-undecided-stimA, use both stimA state (ditch undecided)

%defaults
Fopt.rules = 'strict'; %relaxed/strict
Fopt.Tmax = NaN; %nan, or transition maximum threshold 

if ~isempty(varargin)
    num_args = numel(varargin) / 2;
    fnames = varargin(1:2:end);
    fvals = varargin(2:2:end);
    for idx = 1:num_args
        %check for typo first
        if ~isfield(Fopt,fnames{idx})
            error(sprintf('unknown argument: %s',fnames{idx}))
        end
        Fopt.(fnames{idx}) = fvals{idx};
    end
end

if all(isnan(options.stim_pulse)) %there's a legit pulsing stim 
    sample_duration = options.tmax;
else %constant stim 
    sample_duration = sum(options.stim_pulse);
end

%stimulus_duration = options.stim_pulse(1);
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
%timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});

if isempty(timecourse.state) %just so startsWith() wont error, then exit func
    timecourse.state = '';
end

%just remove stuff before the first valid stay. Must remove
%everything until the first leave following a stay
leave_states = find(strcmp(timecourse.state,'leave'));
first_stay = find(startsWith(timecourse.state,'stim'),1,'first');
if isempty(leave_states)|| isempty(first_stay) || sum(leave_states > first_stay) == 0
    %there's nothing here, abort function call
    Fopt.rules = false;
    %pass empty output
    out_vars = {'event_time','duration','sample_time','sample_number','state','seqID'};
    valid_data = array2table(zeros(0,numel(out_vars)),'VariableNames',out_vars);
else
    
    %finish removing invalid datapoints, run function
    
    %leave states after the first artifical stay
    valid_start = leave_states(leave_states > first_stay);
    valid_start = valid_start(1);
    timecourse(1:valid_start,:) = []; %cut everything before this point
end


switch Fopt.rules
    case 'strict'
        %rule's strict here. Stim A is the reference stimulus, B varies
        %a valid sequence goes [undecided,leave,A,undecided,leave,undecided,B,undecided,leave]
        %where the first, and second stimulus MUST be different. However,
        %can be either A-B, or B-A.
        valid_sequence = {'leave','undecided','stim_','undecided','leave','undecided',...
            'stim_','undecided','leave'}';
        switch options.state_def  %whether simulation aknowledges "undecided states"
            case 'active_states'
                valid_sequence = valid_sequence(~ismember(valid_sequence,'undecided'));
        end
        seq_stim_inds = startsWith(valid_sequence,'stim');
        seq_length = numel(valid_sequence);
        out_vars = {'event_time','duration','sample_time','sample_number','state','seqID'};
        valid_data = array2table(zeros(0,numel(out_vars)),'VariableNames',out_vars);
        Ns = numel(timecourse(:,1)); %total states in the data
        if Ns < seq_length
            stop = []; %skip loop if there's nothing there
        else
            stop = (Ns - seq_length)+1;
        end
        
        %record an index, use to find # A-B vs B-A etc
        seq_counter = array2table(zeros(seq_length,1),'VariableNames',{'seqID'});
        %search as sliding window
        for idx = 1:stop %search as sliding window
            seq_inds = idx:idx+seq_length-1;
            seq = timecourse(seq_inds,:);
            isvalid = sum(cellfun(@startsWith,seq.state,valid_sequence)) == seq_length;
            if isvalid
                %now see if the two stimuli are different
                seq_stims = seq.state(seq_stim_inds);
                stimcheck = numel(unique(seq_stims)) == 2;
                if stimcheck
                    %meets strict rule criteria
                    seq_counter.seqID = seq_counter.seqID + 1; %increment counter
                    seq = [seq,seq_counter]; %cat ID to sequence data
                    valid_data = cat(1,valid_data,seq);
                end
            end
        end
        
        
    case 'relaxed'
        %more chillaxed kinda thing. Stim A is the reference stimulus, B varies
        %a valid sequence goes [A,undecided,leave,undecided,B,undecided,leave]
        %where the first, and second stimulus MUST be different. However,
        %can be either A-B, or B-A. Relaxation piece--- in the sequence
        %leave-A-undecided-A-undecided-leave-B, take duration for A as both
        %A states. Likewise for subsequent B states
        %valid_sequence = {'stim_A','undecided','leave','undecided',...
        %    'stim_B','undecided','leave'};
        leave_states = find(strcmp(timecourse.state,'leave')); %find this again..
        valid_sequence = {'stim_','undecided','leave','undecided',...
            'stim_'}';
        switch options.state_def  %whether simulation aknowledges "undecided states"
            case 'active_states'
                valid_sequence = valid_sequence(~ismember(valid_sequence,'undecided'));
        end
        seq_stim_inds = startsWith(valid_sequence,'stim');
        seq_length = numel(valid_sequence);
        out_vars = {'event_time','duration','sample_time','sample_number','state','seqID'};
        valid_data = array2table(zeros(0,numel(out_vars)),'VariableNames',out_vars);
        Ns = numel(timecourse(:,1)); %total states in the data
        if Ns < seq_length
            stop = []; %skip loop if there's nothing there
        else
            stop = (Ns - seq_length)+1;
        end
        %record an index, use to find # A-B vs B-A etc
        seq_counter = 0; %this will be variable length
        %seq_counter = array2table(zeros(seq_length,1),'VariableNames',{'seqID'}); %this will be variable length
        %search as sliding window
        for idx = 1:stop %search as sliding window
            seq_inds = idx:idx+seq_length-1;
            seq = timecourse(seq_inds,:);
            isvalid = sum(cellfun(@startsWith,seq.state,valid_sequence)) == seq_length;
            if isvalid
                %now see if the two stimuli are different
                seq_stims = seq.state(seq_stim_inds);
                stimcheck = numel(unique(seq_stims)) == 2;
                %ensure there's another subsequent leave state, AND previous one
                leavecheck = sum(leave_states > seq_inds(end)) > 0 && sum(leave_states < idx) > 0;
                if stimcheck && leavecheck
                    %meets alternation & subsequent-leave criteria
                    seq_counter = seq_counter + 1; %increment counter
                    %now we have to search backwards & forwards for more stay-states
                    last_leave = leave_states(leave_states < idx); %timecourse always starts with a leave
                    last_leave = last_leave(end); %... make sure you get the last one dude
                    prev_states = timecourse(last_leave:idx-1,:);
                    %intervening stimuli will always be consistent... dummy
                    next_leave = leave_states(leave_states > seq_inds(end));
                    next_leave = next_leave(1);
                    next_states = timecourse(seq_inds(end)+1:next_leave,:);
                    %add the additional states to sequence
                    seq = [prev_states;seq;next_states];
                    %add in sequence IDs and you're good
                    currIDs = repmat(seq_counter,size(seq,1),1);
                    seq(:,{'seqID'}) = array2table(currIDs);
                    valid_data = cat(1,valid_data,seq);
                end
            end
        end
        
end

Tseq = {'stim_','undecided','leave','undecided','stim_'}'; %transition sequence
switch options.state_def  %whether simulation aknowledges "undecided states"
    case 'active_states'
        Tseq = Tseq(~ismember(Tseq,'undecided'));
        error('you need to fix this part by running data through it...')
end
num_seq = max(valid_data.seqID);
Tvars = {'event_time','u1','leave','u2','total','type'};
Tinfo = array2table(NaN(num_seq,numel(Tvars)),'variablenames',Tvars);
Tinfo.type = cell(size(Tinfo.type));
rec_vars = {'u1','leave','u2'}; %these are the ones where it's just a straight transfer
for idx = 1:num_seq
    seq = valid_data(valid_data.seqID == idx,:);
    leave_states = find(strcmp(seq.state,'leave'));
    if numel(leave_states) ~= 3,error('something wrong with data parsing');end
    Tleave = leave_states(2); %middle one should be the transition 
    Tinfo.event_time(idx) = seq.event_time(Tleave);
    seq = seq(Tleave-2:Tleave+2,:);
    %double check again 
    isvalid = sum(cellfun(@startsWith,seq.state,Tseq)) == numel(Tseq);
    if ~isvalid,error('something wrong with data parsing #2');end
    stim_inds = startsWith(seq.state,'stim_');
    Ttype = seq.state(stim_inds);
    Ttype = strrep(Ttype,'stim_','');
    Ttype = cat(2,Ttype{:});
    Tinfo.type{idx} = Ttype;
    Tinfo.total(idx) = sum(seq.duration(~stim_inds));
    Tinfo{idx,rec_vars} = seq.duration(~stim_inds)';
end

%make sure you didn't get doubled up here somehow 
[~,iU,~] = unique(Tinfo.event_time,'stable');
if numel(iU) ~= num_seq,error('something wrong with data parsing #3');end

if ~isnan(Fopt.Tmax) 
    %impose threshold for stimulus transitions 
    longseq = find(Tinfo.total > Fopt.Tmax); %index does give seqID val 
    longseq = ismember(valid_data.seqID,longseq);
    %remove entries
    valid_data = valid_data(~longseq,:); 
    %do the same for Tinfo 
    Tinfo = Tinfo(Tinfo.total <= Fopt.Tmax,:);
end


%don't include seqID, will always make everything unique
[~,iU,iC] = unique(valid_data(:,1:end-1),'rows','stable');
Nuniq = numel(iU);
if Nuniq ~= size(valid_data,1)
    dups = cellfun(@(x,y) sum(x == iC),num2cell(iC));
    dups = dups > 1; %okay these are the duplicate entries
    %unique always returns the first of any duplicate (with out without 'stable', since valid_data has chrono order)
    %So, the output will have overlapping segments in he correct order, minus the overlap.
    valid_data.seqID = get_overlapIDs(valid_data.seqID,dups); %set all overlapping segments to the same seqID.
    valid_data = valid_data(iU,:); %now cut the overlap
end

%for "decision time"
%dec_time = curr_rec.sample_time(end); %use the sample clock for this, more accurate if there's
% %multiple "samples" within the same pulse i.e. stay-undecided-stay (undecided is usally extremely brief, probably
% %just momentarily below our threshold)
% out.decision_time(idx-1) = dec_time - leave_dur;

end


%code has find_stay_durations() redone with nicer table variables:

%
% switch pref_mode
%     case 'off'
%
%         leave_states = find(strcmp(timecourse.state,'leave'));
%         invalid_switches = zeros(size(leave_states));
%         for idx = 1:numel(leave_states)
%             %leave_durr = timecourse{leave_states(idx),2}; %How long leave-state lasted
%             %leave_start = timecourse{leave_states(idx),1} - leave_durr; %when that leave state started
%             leave_sample = timecourse.sample_number(leave_states(idx));
%             curr_rec = timecourse(1:leave_states(idx),:); %duration record up to that point
%             %check for no prior stay-state (2 sequential leaves after 1st artificial stay)
%             %or for leave states prior to last stay (sequential leaves during normal sim)
%             prev_stay = find(startsWith(curr_rec.state,'stim'), 1, 'last');
%             prev_leave = [0;leave_states]; %just so this wont break with idx == 1
%             prev_leave = max(prev_leave(prev_leave < leave_states(idx)));
%             if isempty(prev_stay) || prev_leave > prev_stay
%                 invalid_switches(idx) = 1;
%             else
%                 prev_stay = curr_rec(prev_stay,:); %stay-state prior to leave-transition
%                 prev_stay_sample = prev_stay.sample_number;
%                 if leave_sample ==  prev_stay_sample
%                     %a stay-state existed in this sample before transitioning to leave
%                     %i.e. legitimate switch in the middle of a sample
%                 else
%                     %it's bogus or something
%                     invalid_switches(idx) = 1;
%                 end
%             end
%         end
%
%         invalid_switches = leave_states(invalid_switches == 1);
%         timecourse(invalid_switches,:) = []; %remove invalid leaves
%         undecided_states = strcmp(timecourse.state,'undecided');
%         timecourse(undecided_states,:) = []; %remove undecided states
%
%     case 'on'
%         %just remove stuff before the first valid stay. Must remove
%         %everything until the first leave following a stay
%         leave_states = find(strcmp(timecourse(:,end),'leave'));
%         first_stay = find(startsWith(timecourse(:,end),'stim'),1,'first');
%         valid_start = leave_states(leave_states > first_stay);
%         valid_start = valid_start(1);
%         timecourse(1:valid_start,:) = []; %cut everything before this point
%        %rule's strict here. Stim A is the reference stimulus, B varies
%         %a valid sequence goes [A,undecided,leave,undecided,B,undecided,leave]
%         valid_sequence = {'stim_A','undecided','leave','undecided',...
%             'stim_B','undecided','leave'};
%         seq_length = numel(valid_sequence);
%         okay_inds = [];
%         keyboard
%
%
%
% end
%
%
% %ok... now for the real results
% leave_states = find(strcmp(timecourse.state,'leave'));
% out = NaN(numel(leave_states),3);
% out = array2table(out,'VariableNames',{'duration','decision_time','samples'});
% switch rule_mode
%     case 'on'
%         %{state.timeidx,state.count,state.sample_clock,record_stim_label});
%         info = cell(numel(leave_states),2); %just return the time index & label for verification
% end
% leave_states = cat(1,0,leave_states); %just for indexing
% for idx = 2:numel(leave_states) %start at 2 b/c you appended that zero for indexing
%     curr_rec = leave_states(idx-1)+1:leave_states(idx);
%     curr_rec = timecourse(curr_rec,:); %stay states up to that point
%     total_duration = curr_rec.duration(1:end-1); %durations of all the stay-states (not inc. leave...)
%     total_duration = sum(total_duration);
%     out.duration(idx-1) = total_duration;
%     %stay_start = curr_rec{end-1,1} - curr_rec{end-1,2}; %when previous stay started
%     %leave_start = curr_rec{end,1} - curr_rec{end,2}; %when this leave state started (aka decision time)
%     leave_dur = curr_rec.duration(end);
%     dec_time = curr_rec.sample_time(end); %use the sample clock for this, more accurate if there's
%     %multiple "samples" within the same pulse i.e. stay-undecided-stay (undecided is usally extremely brief, probably
%     %just momentarily below our threshold)
%     out.decision_time(idx-1) = dec_time - leave_dur;
%     total_samples = numel(unique(curr_rec.sample_number)); %use the unique sample numbers
%     out.samples(idx-1) = total_samples;
%     switch rule_mode
%         case 'on'
%             %previous state's time index & label for verification
%             info(idx-1,:) = table2cell(curr_rec(end-1,{'event_time','state'}));
%     end
% end
%
%
%
% end


