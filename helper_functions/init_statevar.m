function state = init_statevar(celltype,options)
timestep = options.timestep;
%---state tracker-------------
state = struct();
state.stay = logical([1 0]);
state.switch = logical([0 1]);
state.undecided = logical([0 0]);
state.last_leave_end = NaN;
%durations = {}; %record duration time, state/stimulus label
state.now = NaN; %pick one state to start with, add pulse to that pool to be sure (make sure this is consistent)
state.state_def = options.state_def; %whether simulation aknowledges "undecided states" 
state.test_time = options.state_test_time / timestep;
state.test_thresh = options.state_test_thresh;
state.thresh_clock = 0;
state.sample_clock = 0;
state.cut_leave_state = options.cut_leave_state / timestep;
state.stim_labels = {'stim_A','stim_B'};
state.current_stimulus = logical([1 0]); %initialize in stim A
state.count = 0;
switch options.GPU_mdl
    case 'off'
        state.pools2compare = [celltype.pool_stay & celltype.excit,...
            celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
    case 'on'
          state.pools2compare = [celltype.pool_stay,celltype.pool_switch];
          state.pools2compare = state.pools2compare(celltype.excit,:);
end
state.GPU_mdl = options.GPU_mdl;
state.init_check_Lext = options.init_check_Rext * timestep;
state.init_check_stop = options.init_check_tmax / timestep; %minimum time for ready2go check
state.noswitch_timeout = options.noswitch_timeout / timestep;
state.no_dominance_timeout = options.no_dominance_timeout / timestep;
state.no_dom_counter = 0;
state.sample_Estay_offset = options.sample_Estay_offset / timestep; %new: this right here