function ext_spikes = back2stay(ext_spikes,state,celltype)
%if currently in a switch state, blank out all incoming spikes to
%excitatory switch cells. Force network back to a stay state

if all(state.now == state.switch) & state.count >= state.back2stay_min
    ext_spikes(celltype.excit & celltype.pool_switch) = 0;
end