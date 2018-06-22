function [state,durations] = test4switch(Sg,state,durations)

if mean(Sg(state.pools2compare(:,state.now))) * 4 < mean(Sg(state.pools2compare(:,~state.now))) %difference by factor of 4
    %switch!
    if all(state.stay == state.now)%we were just in a stay state
        record_stim_label = state.stim_labels{state.current_stimulus}; %record stimulus label
        state.current_stimulus = ~state.current_stimulus; %switch to the other stimulus for next stay-state
    elseif all(state.switch == state.now) %we were just in a switch state
        record_stim_label = NaN;
    end
    durations{state.now} = vertcat(durations{state.now},{state.count,record_stim_label}); %record duration in num timesteps
    state.now = ~state.now; %switch state
    state.count = 0; %reset the clock
else
    state.count = state.count + 1;
end


% if mean(Sg(state.pools2compare(:,state.now))) * 4 < mean(Sg(state.pools2compare(:,~state.now))) %difference by factor of 4
%     %switch!
%     durations{state.now} = vertcat(durations{state.now},state.count); %record duration in num timesteps
%     state.now = ~state.now; %switch state
%     state.count = 0; %reset the clock
% else
%     state.count = state.count + 1;
% end



