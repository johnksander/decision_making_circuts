function curr_noise = NFSnoise_adjustment(curr_noise,current_info,state,durations,celltype)
%-----noise adjustment for looking at switching dynamics-----
%adjust noise to particular cells in order to force a switch
%50 ms duration, minimum of 1 second into state
%state when noise adjusted depends on the switch being examined
%target cells and noise adjustment depends on particular simulation

%!!!target state hardcoded to stim A!!!

push_noise = false;

%determine if we're in the pre-switch push window
if state.count >= current_info.NFS_onset_min & state.now == state.stay
    %minimum of 1 second into a stay state
    num_stay_states = numel(durations{state.stay});
    if num_stay_states > 2 & mod(num_stay_states,2) == 0
        %last recorded duration was for stim B & skip over first artifically induced A state
        %we're in the stimuli A stay state now, at least 1 second in
        if state.count <= current_info.NFS_stoppush
            %lastly, we havn't been forcing a switch for more than Xms
            push_noise = true;
        end
    end
end


%determine if we're in the post-forced-switch push window
%find out if we've just switched from A to [the switch before B], < Xms ago
if state.count <= current_info.NFS_afterswitch_push & state.now == state.switch
    %if switch Xms ago & we're in the switch state
    num_stay_states = numel(durations{state.stay});
    if num_stay_states > 2 & mod(num_stay_states,2) == 1
        %last recorded duration was for stim A & skip over first artifically induced A state
        %we're in the switch after A now
        last_duration = durations{state.stay};
        last_duration = last_duration(end);
        if last_duration >= current_info.NFS_onset_min & last_duration <= current_info.NFS_stoppush
            %make sure the switch happened during the noise push window
            push_noise = true;
        end
    end
end


if push_noise
    switch current_info.NFS
        case 'Estay'
            targets = celltype.excit & celltype.pool_stay;
            noise_push = -current_info.NFS_noisepush; %decrease stay excitability
        case 'Eswitch'
            targets = celltype.excit & celltype.pool_switch;
            noise_push = current_info.NFS_noisepush; %increase switch excitability
        case 'Istay'
            targets = celltype.inhib & celltype.pool_stay;
            noise_push = current_info.NFS_noisepush; %increase stay inhibition
        case 'Iswitch'
            targets = celltype.inhib & celltype.pool_switch;
            noise_push = -current_info.NFS_noisepush; %decrease switch inhibition
    end
    curr_noise(targets) = curr_noise(targets) + noise_push; %apply noise push
end