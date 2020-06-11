function modelfile = special_ratejob_model(options)
%only used for getting some spikerates retroactively... needed to eliminate 
%durations check to make this happen quickly
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
%seperate for efficiency 
W_Ex = W(:,celltype.excit);
W_In = W(:,celltype.inhib);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell basics---------------
Erev = 0; %reversal potential, excitatory
Irev = -70e-3; %reversal potential, inhibitory
Gg = 10e-9; %max conductance microSiemens
Tsyn_Ex = 50e-3; %excitatory gating time constant, ms
Tsyn_In = 10e-3; %inhibitory gating time constant, ms
El = -70e-3; %leak potential mV
Ek = -80e-3; %potassium potential mV
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%----noisy input---------------
Tau_ext = NaN(pool_options.num_cells,1); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 3.5e-3; % was 2e-3
Tau_ext(celltype.inhib) = 2e-3; % was 5e-3
initGext = 10e-9; %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9; %increase noisy conducrance, nano Siemens
Rext = 1540; %poisson spike train rate for noise, Hz
%----adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----depression----------------
Pr = .1; %release probability
Td.fast = .3;%synaptic depression time constant, seconds (fast)
Td.slow = 10;%slow time constant
fD = .05; %availability factor: (max docked vessicles / max pooled vessicles)
switch options.fastslow_depression
    case 'off'
        fD = 0; Td.slow = Td.fast;
end
%----timecourse----------------
timestep = options.timestep;
num_timepoints = round(options.tmax / timestep) + 1; %same as numel(0:dt:Tmax)
Lext = Rext * timestep; %poisson lambda for noisy conductance

%simulation trial loop
%-------------------------------------------------------------------------
update_logfile(':::Starting simulation:::',options.output_log)
num_trials = numel(options.trial_stimuli(:,1));
sim_results = cell(num_trials,4);

for trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = zeros(pool_options.num_cells,1);
    V = V + El;%inital value of membrane potential is leak potential
    %---stimuli info--------------
    stim_info = init_stimvar(celltype,pool_options,options);
    %for trials, work out passing trial_stimuli indexed by trialidx to init_statevar()
    %---noisy conductance---------
    ext_inds.I = 1;
    ext_inds.E = 2;
    Gext_Ex = initGext + zeros(pool_options.num_cells,1); %initialize at leak conductance
    Gext_In = initGext + zeros(pool_options.num_cells,1);
    %---adaptation conductance----
    Gsra = zeros(size(V));
    %---gating & depression-------
    Sg_Ex = zeros(sum(celltype.excit),1);
    Sg_In = zeros(sum(celltype.inhib),1);
    Dfast = ones(size(V)); %synaptic depression: fast
    Dslow = ones(size(V)); %synaptic depression: slow
    %---spikes--------------------
    switch options.record_spiking
        case 'on'
            spikes = zeros(pool_options.num_cells,num_timepoints); %preallocating the whole thing in this one...
        otherwise %see if a smaller matrix needs to be allocated for ratelim check
            switch options.ratelim.check
                case 'on'
                    RLCstart = options.ratelim.start /timestep; RLCstop = options.ratelim.stop /timestep;
                    spikes = zeros(pool_options.num_cells,RLCstop-RLCstart);
                    options.record_spiking = 'ratelim_only'; %this will reset -> 'off' after check
            end
    end
    %---state tracker-------------
    durations = {}; %record duration time, state/stimulus label
    state = init_statevar(celltype,options);
    state.now = state.undecided; %this will always be true when V init to El
    %---last init-----------------
    experiment_set2go = false; %when experiment is ready to go
    avail_noise.Estay = 1; avail_noise.Eswitch = 1;
    timepoint_counter = 1;
    
    while timepoint_counter < num_timepoints
        
        timepoint_counter = timepoint_counter+1;
        state.timeidx = timepoint_counter; %just so I don't have to pass a million things...
        
        WS_Ex = W_Ex*Sg_Ex;
        WS_In = W_In*Sg_In;
        
        E_diff = Erev-V;
        I_diff = Irev-V;
        I = E_diff.*WS_Ex.*Gg + I_diff.*WS_In.*Gg + Gext_Ex.*E_diff + Gext_In.*I_diff;
        dVdt = ((El - V + (delta_th.*exp((V - Vth)./delta_th)))./Rm) + (Gsra.*(Ek - V)) + I;
        V = ((dVdt./Cm) .* timestep) + V;
                 
        Gext_Ex = Gext_Ex - (Gext_Ex./Tau_ext) .* timestep;%noisy conductance
        Gext_In = Gext_In - (Gext_In./Tau_ext) .* timestep;
        Gsra = Gsra - (Gsra./Tsra) .* timestep;%adaptation conductance
        Sg_Ex = Sg_Ex - (Sg_Ex./Tsyn_Ex) .* timestep;%synaptic gating
        Sg_In = Sg_In - (Sg_In./Tsyn_In) .* timestep;
        Ddiff = Dslow - Dfast;
        Dfast = Dfast + ((Ddiff ./ Td.fast) .* timestep);%fast syn. depression
        Dslow = Dslow + timestep .* ( ((1 - Dslow) ./ Td.slow) ...
            - fD.* (Ddiff ./ Td.fast)  ); %slow vessicle replacement
       
        spiking_cells = V > spike_thresh;
        if sum(spiking_cells) > 0
            Gsra(spiking_cells) = Gsra(spiking_cells) + detlaGsra; %adaptation conductance
            %synaptic gating, updates with Dfast vessicles
            Espike = spiking_cells(celltype.excit);
            Ispike = spiking_cells(celltype.inhib);
            Sg_Ex(Espike) = Sg_Ex(Espike) + (Pr.* Dfast(spiking_cells & celltype.excit).*(1-Sg_Ex(Espike)));
            Sg_In(Ispike) = Sg_In(Ispike) + (Pr.* Dfast(spiking_cells & celltype.inhib).*(1-Sg_In(Ispike)));
            %depression update (docked vessicles released)
            Dfast(spiking_cells) = Dfast(spiking_cells) - (Pr.* Dfast(spiking_cells));
            V(spiking_cells) = Vreset;
            switch options.record_spiking
                case 'on'
                    spikes(spiking_cells,timepoint_counter) = 1;
                case 'ratelim_only' %this gets set to 'off' after the check
                    if timepoint_counter >=  RLCstart
                        spikes(spiking_cells,timepoint_counter-RLCstart+1) = 1;
                    end
            end
        end
        
        %test for state transition & determine stim availability
        if experiment_set2go
            [state,durations] = test4switch(Sg_Ex,state,durations);
            [state,avail_noise] = check_noise_avail(stim_info,state);
        else
            state.count = state.count + 1; %if you don't run test4switch(), must update this counter outside
        end
        
        %run the bistability check
        if ~experiment_set2go %during bistability check, check_bistability() handles pulse input spikes
            [BScheck,Pspikes,state] = check_bistability(Sg_Ex,state);
            %add pulse spikes (same as below)
            Gext_Ex(celltype.excit) = Gext_Ex(celltype.excit) + deltaGext.*Pspikes;
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
        ext_spikes = poissrnd(Lext,pool_options.num_cells,2);
        %adjust noise input if needed
        if avail_noise.Estay ~= 1
            ext_spikes(celltype.excit & celltype.pool_stay,:) = ...
                poissrnd(avail_noise.Estay.*Lext,sum(celltype.excit & celltype.pool_stay),2);
        end
        if avail_noise.Eswitch ~= 1
            ext_spikes(celltype.excit & celltype.pool_switch,:) = ...
                poissrnd(avail_noise.Eswitch.*Lext,sum(celltype.excit & celltype.pool_switch),2);
        end
        
        %---spiking input from stimulus
        if experiment_set2go
            stim_spikes = timepoint_stimulus(stim_info,state); %get stimulus spikes
            %always exitatory, add 'em both together for one calculation
            ext_spikes(:,ext_inds.E) = ext_spikes(:,ext_inds.E) + stim_spikes;
        end
        
        %update Gexternal. Don't have to index, they get an increase or zero        
        Gext_Ex = Gext_Ex + deltaGext.*ext_spikes(:,ext_inds.E);
        Gext_In = Gext_In + deltaGext.*ext_spikes(:,ext_inds.I);
        
        switch options.fastslow_depression
            case 'off'
                Dslow(:,idx) =  1; %keep constant at 1.
        end
        
        %check timeout for non-switching or non-dominance
        if state.count >= state.noswitch_timeout || state.no_dom_counter >= state.no_dominance_timeout
            update_logfile(':::Bistability check failure:::',options.output_log)
            TOF = timepoint_counter*timestep;
            update_logfile(sprintf('---no switch/dominance timeout at t=%.2f(s)',TOF),options.output_log)
            return
        end
        
        %rate limit check
        switch options.ratelim.check
            case 'on'
                if timepoint_counter == RLCstop
                    [options,Rcheck] = check_rate_limit(spikes,celltype,options);
                    switch Rcheck.status
                        case 'fail'
                            return
                        case 'pass'
                            update_logfile('---passed rate limit check',options.output_log)
                    end
                    switch options.record_spiking
                        case 'off' %clear spike matrix from memory
                            clear spikes
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
    %if isempty(durations),return;end
    %trim_Bcheck = find(startsWith(durations(:,end),'leave'), 1, 'first');
    %durations = durations(trim_Bcheck+1:end,:); if isempty(durations),return;end
    sim_results{trialidx,1} = durations;
    
    %record ratelim check's rough spikerate estimate
    switch options.ratelim.check
        case 'on'
            sim_results{trialidx,4} = Rcheck;
    end
    
    switch options.record_spiking
        case 'on'
            
            %---switch data recording-----
            [~,valid_events] = find_stay_durations(durations,options,'verify');
            good_data = cat(1,valid_events.spiking_data{:});
            valid_events = valid_events.event_time(good_data);
            valid_events = cat(1,valid_events{:}); %event times for events with good spiking data
            all_events = cat(1,durations{:,1}) .* options.timestep;
            valid_events = ismember(all_events,valid_events);
            switch_record = durations(valid_events,:);
            
            %check if no switches met recording criteria, terminate if needed
            if isempty(switch_record) | numel(switch_record(:,1)) < 1,return;end
            
            %now get these spike timecourses and save them
            num_switches = numel(switch_record(:,1)); %cannot believe this var name is still free
            num_preswitch_samples = round(options.record_preswitch/timestep); %window samples
            num_postswitch_samples = round(options.record_postswitch/timestep);
            num_switch_samples = num_preswitch_samples + num_postswitch_samples;
            spiking_output = NaN(pool_options.num_cells,num_switch_samples,num_switches);
            for tc_idx = 1:num_switches
                record_win = switch_record{tc_idx,1};
                record_win = record_win-(num_preswitch_samples-1):record_win+num_postswitch_samples;
                spiking_output(:,:,tc_idx) = spikes(:,record_win);
            end
            
            sim_results{trialidx,2} = spiking_output;
            sim_results{trialidx,3} = switch_record; %save this thing along with it
    end
end

update_logfile('---Simulation complete---',options.output_log)
savename = fullfile(options.save_dir,options.sim_name);
save(savename,'sim_results','options')