function state = init_statevar(celltype,options)
timestep = options.timestep;
%---state tracker-------------
state = struct();
state.stay = logical([1 0]);
state.switch = logical([0 1]);
state.undecided = logical([0 0]);
state.last_leave_end = NaN;
durations = {}; %record duration time, state/stimulus label
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
state.pools2compare = [celltype.pool_stay & celltype.excit,...
    celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
state.init_check_Lext = options.init_check_Rext * timestep;
state.init_check_stop = options.init_check_tmax / timestep; %minimum time for ready2go check
state.noswitch_timeout = options.noswitch_timeout / timestep;
state.sample_Estay_offset = options.sample_Estay_offset / timestep; %new: this right here