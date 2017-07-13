function stim_spikes = timepoint_stimulus(stim_info,state)
%inputs: stimuli info & state parameter structures

stim_spikes = zeros(stim_info.num_cells,1);

if sum(state.switch == state.now) == numel(state.now)
    %we're in a switch state, don't do anything
elseif sum(state.stay == state.now) == numel(state.now) %we're in a stay state
    %we're in a stay state, do something
    %add stimulus-specific spiketrain 
    if strcmp(state.stim_labels{state.current_stimulus},'A')
        stim_spikes(stim_info.targ_cells) =  poissrnd(stim_info.stimA_lambda,sum(stim_info.targ_cells),1);
    else %strcmp(state.stim_labels{state.current_stimulus},'B') 
        stim_spikes(stim_info.targ_cells) =  poissrnd(stim_info.stimB_lambda,sum(stim_info.targ_cells),1);        
    end
end






