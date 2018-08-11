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
            t = mod(state.count,seq_length); %find if we're in the on, or off part of that sequence
            if t <= stim_info.pulse(1)
                %we're in the initial "stimulus on" part
                apply_stimulus = true;
            else
                %we're past the "stimulus on" part, the stimulus is off
                apply_stimulus = false;
            end
    end
else
    %we're in a switch or undecided state, don't do anything
    %elseif all(state.switch == state.now) || all(state.undecided == state.now)
    apply_stimulus = false;
end



if apply_stimulus
    %add stimulus-specific spiketrain
    if strcmp(state.stim_labels{state.current_stimulus},'A')
        stim_spikes(stim_info.targ_cells) =  poissrnd(stim_info.stimA_lambda,sum(stim_info.targ_cells),1);
    else %strcmp(state.stim_labels{state.current_stimulus},'B')
        stim_spikes(stim_info.targ_cells) =  poissrnd(stim_info.stimB_lambda,sum(stim_info.targ_cells),1);
    end
end



