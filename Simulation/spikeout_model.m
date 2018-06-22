function modelfile = spikeout_model(options)
%this model func designed to save whole spiking matrix, for short jobs 
modelfile = mfilename; %for backup purposes


%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
%--------------------------------------------------------------------------
%build the connectivity matrix
%make celltype logicals
celltype = celltype_logicals(pool_options);
%make these vectors into logical matricies
[EtoE,EtoI,ItoE] = connection_logicals(celltype,pool_options.num_cells);
%make connection scheme based off connection matricies
connection_scheme = EtoE | EtoI | ItoE;  %plan is:  EtoE | EtoI | ItoE;
%make synaptic connection matrix
W = rand(pool_options.num_cells) < pool_options.p_conn; %connection matrix
W = double(W & connection_scheme); %filter out connections not allowed by scheme
%modify connection weights
W(W > 0 & EtoE) = options.EtoE;
W(W > 0 & ItoE) = options.ItoE;
W(W > 0 & EtoI) = options.EtoI;
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential, excitatory
Erev(celltype.inhib) = -70e-3; %reversal potential, inhibitory
Gg = NaN(pool_options.num_cells,1); %max connductance
Gg(celltype.excit) = 10e-9; %excitatory max connductance microSiemens
Gg(celltype.inhib) = 10e-9; %inhibitory max connductance microSiemens
Pr = NaN(pool_options.num_cells,1); %release probability
Pr(celltype.excit) = .2; %excitatory release probability
Pr(celltype.inhib) = .2; %inhibitory release probability
Td = .3;%synaptic depression time constant, seconds
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50e-3; %excitatory gating time constant, ms
Tsyn(celltype.inhib) = 10e-3; %inhibitory gating time constant, ms
%----cell basics--------------
El = -70e-3; %leak potential mV
Ek = -80e-3; %potassium potential mV
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%----noisy input--------------
Tau_ext = NaN(pool_options.num_cells,1); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 2e-3;
Tau_ext(celltype.inhib) = 5e-3;
initGext = 10e-9; %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9; %increase noisy conducrance, nano Siemens
Rext = 1400; %poisson spike train rate for noise, Hz
%---adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----timecourse---------------
timestep = .25e-3; %.25 milisecond timestep
timevec = 0:timestep:options.tmax;
num_timepoints = numel(timevec);
Lext = Rext * timestep; %poisson lambda for noisy conductance

%simulation trial loop
%-------------------------------------------------------------------------
update_logfile(':::Starting simulation:::',options.output_log)
num_trials = numel(options.trial_stimuli(:,1));
sim_results = cell(num_trials,3);

for trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = NaN(pool_options.num_cells,2);
    V(:,1) = El; %inital value of membrane potential is leak potential
    %---stimuli info--------------
    stim_info = struct();
    switch options.stim_targs
        case 'Eswitch'
            stim_info.targ_cells = celltype.excit & celltype.pool_switch; %Eswitch cells
        case 'Estay'
            stim_info.targ_cells = celltype.pool_stay & celltype.excit; %Estay
        case 'baseline'
            stim_info.targ_cells = logical(zeros(pool_options.num_cells,1)); %no targets
    end
    stimA = options.trial_stimuli(trialidx,1);
    stimB = options.trial_stimuli(trialidx,2);
    stim_info.stimA_lambda = stimA * timestep; %poisson lambda for stimulus conductance
    stim_info.stimB_lambda = stimB * timestep;
    stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
    if all(~isnan(options.stim_pulse))
        stim_info.delivery = 'pulse';
        stim_info.pulse = options.stim_pulse ./ timestep;
    else 
        stim_info.delivery = 'constant';
    end
    %---noisy conductance---------
    Gext = NaN(size(V)); %noisy conductance
    Gext(:,1) = initGext; %initialize at leak conductance
    %---adaptation conductance----
    Gsra = NaN(size(V));
    Gsra(:,1) = 0;
    %---gating & depression-------
    Sg = NaN(size(V)); %synaptic gating
    Sg(:,1) = 0; %initalize at zero??
    D = NaN(size(V)); %synaptic depression
    D(:,1) = 1; %initalize at one
    %---spikes--------------------
    spikes = zeros(pool_options.num_cells,num_timepoints); %preallocating the whole thing in this one... 
    %---state tracker-------------
    state = struct();
    state.stay = logical([1 0]);
    state.switch = logical([0 1]);
    durations = cell(1,2); %record duration times (stay, switch)
    state.now = state.stay; %pick one state to start with, add pulse to that pool to be sure (make sure this is consistent)
    state.stim_labels = {'A','B'};
    state.current_stimulus = logical([1 0]); %initialize in stim A
    state.count = 0;
    state.pools2compare = [celltype.pool_stay & celltype.excit,...
        celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
    %state.ready_mintime = .4 / timestep; %minimum time for ready2go check
    state.init_check_Lext = options.init_check_Rext * timestep;
    state.init_check_stop = options.init_check_tmax / timestep; %minimum time for ready2go check
    state.noswitch_timeout = options.noswitch_timeout / timestep;
    if options.force_back2stay
        state.back2stay_min = .5 / timestep; %wait this long before forcing switch back
    end
    %---switch data recording-----
    %250ms before switch, 150ms after
    num_preswitch_samples = 250e-3/timestep;
    num_postswitch_samples = 150e-3/timestep;
    num_switch_samples = num_preswitch_samples + num_postswitch_samples;
    switch_record = {};
    %---last init-----------------
    experiment_set2go = false; %when experiment is ready to go
    undecided_state = false; %in between stimulus delivery pulses in stay state
    timepoint_counter = 1;
    idx = 2; %keep indexing vars with idx fixed at 2
    
    while timepoint_counter <= num_timepoints
        
        timepoint_counter = timepoint_counter+1;
        
        %get current for this timepoint
        stim_info.timeidx = timepoint_counter; %just so I don't have to pass a million things...
        
        %loop equations
        I = (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
        I = I + (unique(Erev(celltype.inhib)) - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*unique(Gg(celltype.inhib));
        dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm)...
            + (Gsra(:,idx-1).*(Ek-V(:,idx-1)))...
            + (Gext(:,idx-1).*(Erev-V(:,idx-1))) + I;
        
        V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
        
        Gsra(:,idx) = Gsra(:,idx-1) - ((Gsra(:,idx-1)./Tsra) .* timestep); %adaptation conductance
        Gext(:,idx) = Gext(:,idx-1) - ((Gext(:,idx-1)./Tau_ext) .* timestep); %noisy conductance
        D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
        Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
        
        spiking_cells = V(:,idx) > spike_thresh;
        if sum(spiking_cells) > 0
            spiking_cells = V(:,idx) > spike_thresh;
            Gsra(spiking_cells,idx) = Gsra(spiking_cells,idx) + detlaGsra; %adaptation conductance
            Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
                (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
            D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
            V(spiking_cells,idx) = Vreset;
            spikes(spiking_cells,timepoint_counter) = 1;
        end
        
        %make sure we're not in an undecided state (pulse stimulus delivery)
        if experiment_set2go
            undecided_state = check_undecided_state(stim_info,state);
        end
        
        %test for state transition
        if ~undecided_state
            [state,durations] = test4switch(Sg(:,idx),state,durations);
        else
            state.count = state.count + 1; %if you don't run test4switch(), must update this counter outside 
        end
        
        %run the bistability check
        if ~experiment_set2go %during bistability check, check_bistability() handles pulse input spikes
            [BScheck,Pspikes] = check_bistability(Sg(:,idx),state,durations);
            %add pulse spikes (same as below)
            Gext(:,idx) = Gext(:,idx) + (deltaGext.*Pspikes);
            switch BScheck.status
                case 'fail'
                    update_logfile(':::Bistability check failure:::',options.output_log)
                    update_logfile(sprintf('---at t=%.2f(s)',timepoint_counter*timestep),options.output_log)
                    update_logfile(BScheck.Fcode,options.output_log)
                    return
                case 'pass'
                    update_logfile('---passed bistability check',options.output_log)
                    experiment_set2go = true; %we're ready to roll
            end
        end
        
        %input spikes: noise & stimulus
        %---noisy spiking input from elsewhere
        ext_spikes = poissrnd(Lext,pool_options.num_cells,1);
        if undecided_state
             %we're in between stimulus delivery pulses
             %do half-noise to all E-cells-- gives barely spiking undecided state
             ext_spikes(celltype.excit) = poissrnd(Lext*.5,sum(celltype.excit),1); %half noise E-cells
        end
        %---spiking input from stimulus  
        if experiment_set2go 
            stim_spikes = timepoint_stimulus(stim_info,state); %get stimulus spikes
            ext_spikes = ext_spikes + stim_spikes; %add 'em both together for one calculation
        end
        %if in switch state, force back to a stay state (after back2stay_min time)
        if options.force_back2stay
            ext_spikes = back2stay(ext_spikes,state,celltype); 
        end
        %update Gexternal. Don't have to index, they get an increase or zero
        Gext(:,idx) = Gext(:,idx) + (deltaGext.*ext_spikes); 
        
        %lag equation vars for next timepoint
        V = next_timepoint(V);
        Gsra = next_timepoint(Gsra);
        Gext = next_timepoint(Gext);
        D = next_timepoint(D);
        Sg = next_timepoint(Sg);
        
        %check timeout for non-switching
        if state.count >= state.noswitch_timeout
            update_logfile(':::Bistability check failure:::',options.output_log)
            TOF = timepoint_counter*timestep;
            update_logfile(sprintf('---no switch timeout at t=%.2f(s)',TOF),options.output_log)
            return
        end
        
        %check if we need to record a switch index
        if experiment_set2go
            %find out if we've just switched from A to [the switch before
            %B] Xms ago (or vice-verse, this is agnostic to specific stimuli)
            if state.count == num_postswitch_samples && all(state.now == state.switch)
                %if switch Xms ago & we're in the switch state
                prev_state = durations{state.stay}; %find info on the stay state we just switched out of
                %get duration, make sure it's > num_preswitch_samples (don't want to see another switch in there)
                last_duration = prev_state{end,1};
                switch stim_info.delivery %or see undecided state (non)activity, in pulse scenario
                    case 'pulse' %see how long you've actually been in decided state (mod(t,squence time))
                      last_duration = mod(last_duration,sum(stim_info.pulse));
                end
                if last_duration > num_preswitch_samples
                    %record the switch time, stimulus type, state duration
                    prev_ST = timepoint_counter - state.count;
                    prev_stim = prev_state{end,2}; %find the last stimulus type
                    switch_record = vertcat(switch_record,{prev_ST,prev_stim,last_duration}); 
                end
            end
        end
        
        %progress tracking...
        if mod(timepoint_counter,floor(num_timepoints * .05)) == 0 %5 percent
            progress = (timepoint_counter /  num_timepoints) * 100;
            message = sprintf('Simulation %.1f percent complete',progress);
            update_logfile(message,options.output_log)
        end
    end
    
    %remove the first artificially induced stay state & subsequent switch state
    durations = cellfun(@(x) x(2:end,:),durations,'UniformOutput',false);
    sim_results{trialidx,1} = durations;
    
    %also remove the first artificially induced stay state from this record
    switch_record = switch_record(2:end,:);
    %now get these spike timecourses and save them 
    num_switches = numel(switch_record(:,1)); %cannot believe this var name is still free
    spiking_output = NaN(pool_options.num_cells,num_switch_samples,num_switches);
    for tc_idx = 1:num_switches
        record_win = switch_record{tc_idx,1};
        record_win = record_win-(num_preswitch_samples-1):record_win+num_postswitch_samples;
        spiking_output(:,:,tc_idx) = spikes(:,record_win);
    end
    sim_results{trialidx,2} = spiking_output;
    sim_results{trialidx,3} = switch_record; %save this thing along with it
end
update_logfile('---Simulation complete---',options.output_log)
savename = fullfile(options.save_dir,options.sim_name);
save(savename,'sim_results','options')


