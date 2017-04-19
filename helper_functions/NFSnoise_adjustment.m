function curr_noise = NFSnoise_adjustment(curr_noise,current_info,state,durations,celltype)
%-----noise adjustment for looking at switching dynamics-----
%adjust noise to particular cells in order to force a switch 
%50 ms duration, minimum of 1 second into state
%state when noise adjusted depends on the switch being examined
%target cells and noise adjustment depends on particular simulation 

%!!!target state hardcoded to stim A!!!

if state.count >= current_info.NFS_onset_min & state.now == state.stay
    %minimum of 1 second into a stay state
    num_stay_states = numel(durations{state.stay});
    if num_stay_states > 2 & mod(num_stay_states,2) == 0
        %last recorded duration was for stim B & skip over first artifically induced A state
        %we're in the stimuli A stay state now, at least 1 second in
        if state.count <= current_info.NFS_stoppush
            %lastly, we havn't been forcing a switch for more than 50ms 
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
    end
end