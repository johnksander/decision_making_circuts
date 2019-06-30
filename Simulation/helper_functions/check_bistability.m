function [BScheck,Pspikes,state] = check_bistability(Sg,state)

BScheck.status = 'running'; %initalize
target_cells = state.pools2compare(:,state.stay);
Pspikes = zeros(size(target_cells));
%see if we're still in the bistability check window
check_window = state.count <= state.init_check_stop;
%add another one for "check over", just in case isnan(state.count) or something...
check_over = state.count >= state.init_check_stop;
%give pulse up until 2/3 into check window, then stop
pulse_window = state.count <= floor(state.init_check_stop *8); % *.66 would give 2/3
%if the pulse window is 1/3 over, start checking state dominance
check_dom = state.count >= floor(state.init_check_stop *.33) && check_window;


if pulse_window
    %give pulse spikes
    Pspikes(target_cells) = ...
        poissrnd(state.init_check_Lext,sum(target_cells),1); %external spikes to noise conductance
end

if check_dom
    %test if there's an active state, check if it's the targeted cells (E-stay)
    Smu = [Sg(state.pools2compare(:,1)), Sg(state.pools2compare(:,2))];
    Smu = mean(Smu);
    active_state = Smu - fliplr(Smu); %A-B, B-A
    state.now = active_state > state.test_thresh; %will return undecided, if neither active
    
    targ_active = all(state.now == state.stay); %if it's the state we intended
    if ~targ_active
        state.thresh_clock = state.thresh_clock + 1; 
    elseif targ_active
        state.thresh_clock = 0; %reset the clock
    end
end

%network shouldn't leave target state during the check period
no_switches = state.thresh_clock < state.test_time;


%---outcomes---

if check_dom && ~no_switches
    %active network state isn't the targeted cells---- fail
    if ~strcmp(BScheck.status,'fail')
        %make sure hasn't already failed, give correct message
        BScheck.status = 'fail';
        muA = mean(Sg(state.pools2compare(:,state.stay)));
        muB = mean(Sg(state.pools2compare(:,state.switch)));
        BScheck.Fcode = sprintf('---non-dominance (muStay=%.2f,muLeave=%.2f)',muA,muB);
    end
end

%this was an outcome under the old switch-testing scheme 
% if ~no_switches && check_window
%     %active network state isn't the targeted cells during the pulse window---- fail
%     BScheck.status = 'fail';
%     BScheck.Fcode = '---N > 1 switches';
% end

no_failure = ~strcmp(BScheck.status,'fail');

if check_over && no_switches && no_failure
    %should be good, don't forget to renew your sticker next year
    BScheck.status = 'pass';
end






