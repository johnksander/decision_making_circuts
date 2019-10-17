function [state,durations] = test4switch(Sg,state,durations)

Smu = [Sg(state.pools2compare(:,1)), Sg(state.pools2compare(:,2))];
Smu = sum(Smu) ./ numel(Smu(:,1)); %faster than mean()
active_state = Smu - flip(Smu,2); %A-B, B-A
active_state = active_state > state.test_thresh; %will return undecided, if neither active

%non-dominance check
if all(active_state) || all(~active_state) %both or neither 
    state.no_dom_counter = state.no_dom_counter + 1;
else 
    state.no_dom_counter = 0; %reset the clock 
end


switch state.state_def  %whether simulation aknowledges "undecided states" 
    case 'active_states'
    %must be a different state, that's not undecided
    state_transition = any(active_state ~= state.now) && any(active_state);
    case 'include_undecided'
    state_transition = any(active_state ~= state.now);   
end


%check if there's suprathreshold change in state activity 
if state_transition
    state.thresh_clock = state.thresh_clock + 1;
else 
    state.thresh_clock = 0; %reset the clock 
end

%over Xms suprathreshold: switch! (or it's undecided)
if state.thresh_clock == state.test_time 
    
    if all(state.stay == state.now)%we were just in a stay state
        record_stim_label = state.stim_labels{state.current_stimulus}; %record stimulus label
        
    elseif all(state.switch == state.now) %we were just in a switch state
        record_stim_label = 'leave';
        state.last_leave_end = state.timeidx;
        %note--- this will need to be examined more closely if you start having more than one stim again
        state.current_stimulus = ~state.current_stimulus; %switch to the other stimulus for next stay-state
        
    elseif all(state.undecided == state.now) %we were just undecided
        record_stim_label = 'undecided';
        
    end
    durations = vertcat(durations,...
        {state.timeidx,state.count,state.sample_clock,record_stim_label});
    state.now = active_state; %specify new state
    state.count = 0; %reset the clock
    state.thresh_clock = 0; 

else
    state.count = state.count + 1;
end

%under the old scheme 
% 
% if mean(Sg(state.pools2compare(:,state.now))) * 4 < mean(Sg(state.pools2compare(:,~state.now))) %difference by factor of 4
%     %switch!
%     if all(state.stay == state.now)%we were just in a stay state
%         record_stim_label = state.stim_labels{state.current_stimulus}; %record stimulus label
%         state.current_stimulus = ~state.current_stimulus; %switch to the other stimulus for next stay-state
%     elseif all(state.switch == state.now) %we were just in a switch state
%         record_stim_label = 'leave';
%     end
%     durations{state.now} = vertcat(durations{state.now},{state.count,record_stim_label}); %record duration in num timesteps
%     state.now = ~state.now; %switch state
%     state.count = 0; %reset the clock
% else
%     state.count = state.count + 1;
% end
% 

