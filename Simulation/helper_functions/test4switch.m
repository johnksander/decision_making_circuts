function [state,durations] = test4switch(Sg,state,durations)

if mean(Sg(state.pools2compare(:,state.now))) * 4 < mean(Sg(state.pools2compare(:,~state.now))) %difference by factor of 4
    %switch!
    durations{state.now} = vertcat(durations{state.now},state.count); %record duration in num timesteps
    state.now = ~state.now; %switch state
    state.count = 0; %reset the clock
else
    state.count = state.count + 1;
end



