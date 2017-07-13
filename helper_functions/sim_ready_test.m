function [state,durations,experiment_set2go] = sim_ready_test(Sg,state,durations,experiment_set2go)
%see if the simuluation is in the correct state to begin the experiment
%DEPRECIATED!!!!

%do the normal switching business
if mean(Sg(state.pools2compare(:,state.now))) * 4 < mean(Sg(state.pools2compare(:,~state.now))) %difference by factor of 4
    %switch!
    if sum(state.stay == state.now) == numel(state.now) %we were just in a stay state
        record_stim_label = state.stim_labels{state.current_stimulus}; %record stimulus label
        state.current_stimulus = ~state.current_stimulus; %switch to the other stimulus for next stay-state
    else %we were just in a switch state
        record_stim_label = NaN;
    end
    durations{state.now} = vertcat(durations{state.now},{state.count,record_stim_label}); %record duration in num timesteps
    state.now = ~state.now; %switch state
    state.count = 0; %reset the clock
else
    state.count = state.count + 1;
end


%test and see if we're ready to begin the experiment
if sum(state.stay == state.now) == numel(state.now) %we're currenty in the stay state
    if strcmp(state.stim_labels{state.current_stimulus},'B') && state.count  >= state.ready_mintime
        %if we're doing stimulus B right now, and we've been doing it for a little while
            durations = cell(1,2); %blank out the recorded durations
            experiment_set2go = true; %we're ready to roll
    end
end



