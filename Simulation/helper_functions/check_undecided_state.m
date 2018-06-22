function outcome = check_undecided_state(stim_info,state)
%return true if the stimulus delivery method is "pulse"
%and we're currently in the "stimulus off" part of that delivery squence.
%(must be in the stay-state as well, for clarity). 
%otherwise return false for any other scenario 

outcome = false;

switch stim_info.delivery
    case 'pulse'
        if all(state.stay == state.now) 
            %we're in a stay state
            seq_length = sum(stim_info.pulse); %how long for a single on, off sequence
            t = mod(state.count,seq_length); %find if we're in the on, or off part of that sequence
            if t > stim_info.pulse(1)
                %we're past the "stimulus on" part, the stimulus is off
                outcome = true;
            end
        end
end