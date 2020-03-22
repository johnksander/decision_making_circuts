function stim_spikes = timepoint_stimulus(stim_info,state)
%inputs: stimuli info & state parameter structures

stim_spikes = zeros(stim_info.num_cells,1);

%figure out if we need to apply a stimulus
if all(state.stay == state.now) %we're in a stay state
    %we're in a stay state, do something
    switch stim_info.delivery
        case 'constant'
            apply_stimulus = true;
        case 'pulse'
            seq_length = sum(stim_info.pulse); %how long for a single on, off sequence
            t = mod(state.sample_clock,seq_length); %find if we're in the on, or off part of that sequence
            if t <= stim_info.pulse(1) && state.sample_clock > 0
                %we're in the initial "stimulus on" part
                apply_stimulus = true;
            else
                %we're past the "stimulus on" part, the stimulus is off. 
                %Or the sample clock is negative, we're in the delay after leave 
                apply_stimulus = false;
            end
    end
else
    %we're in a switch or undecided state, don't do anything
    %elseif all(state.switch == state.now) || all(state.undecided == state.now)
    apply_stimulus = false;
end



if apply_stimulus
    
    curr_stim = state.stim_labels{state.current_stimulus};
    lambda = stim_info.(curr_stim); %1 x Ntarg vector 
    
    for idx = 1:stim_info.num_targs    
        targs = stim_info.targ_cells(:,idx); %Ncell x Ntarg matrix 
        %add stimulus-specific spiketrain
        stim_spikes(targs) = poissrnd(lambda(idx),sum(targs),1);     
    end
end



