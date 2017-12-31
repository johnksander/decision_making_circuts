function [BScheck,Pspikes] = check_bistability(Sg,state,durations)

BScheck.status = 'running'; %initalize
target_cells = state.pools2compare(:,state.now);
Pspikes = zeros(size(target_cells));

%see if we're still in the bistability check window
check_window = state.count <= state.init_check_stop;
%add another one for "check over", just in case isnan(state.count) or something...
check_over = state.count >= state.init_check_stop;
%give pulse up until 2/3 into check window, then stop
pulse_window = state.count <= floor(state.init_check_stop *1); % *66 would give 2/3
%if the pulse window is 1/3 over, start checking state dominance
check_dom = state.count >= floor(state.init_check_stop *.33) && check_window;
%should have been zero switches during the check period
no_switches = numel(vertcat(durations{:})) == 0;


if pulse_window
    %give pulse spikes
    Pspikes(target_cells) = ...
        poissrnd(state.init_check_Lext,sum(target_cells),1); %external spikes to noise conductance
end

if check_dom
    %see if the current state is actually on & more active than the other
    targ_active = mean(Sg(state.pools2compare(:,state.now))) > 4 * mean(Sg(state.pools2compare(:,~state.now)));
    %test4switch() just checks if the opposite is true, sort of 
end

%---outcomes---

if ~no_switches && check_window
    %states have flipped during the pulse window---- fail
    BScheck.status = 'fail';
    BScheck.Fcode = '---N > 1 switches';
end

if check_dom && ~targ_active
    %targeted cells aren't actually active, or their activity isn't much higher
    %than cells for the other state---- fail
    if ~strcmp(BScheck.status,'fail')
        %make sure hasn't already failed, give correct message
        BScheck.status = 'fail';
        muA = mean(Sg(state.pools2compare(:,state.now)));
        muB = mean(Sg(state.pools2compare(:,~state.now)));
        BScheck.Fcode = sprintf('---non-dominance (muA=%.2f,muB=%.2f)',muA,muB);
    end
end

no_failure = ~strcmp(BScheck.status,'fail');

if check_over && no_switches && no_failure
    %should be good, don't forget to renew your sticker next year
    BScheck.status = 'pass';
end






