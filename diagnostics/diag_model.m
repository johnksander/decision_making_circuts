function [MF,oflag] = diag_model(options)
%this model func designed to save whole spiking matrix, for short jobs
MF = mfilename; %for backup purposes

oflag = false;
%enable recording for... D,Vm,Spikes,(Sg?)  spikes was already good-to-go

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
    
    %need to preallocate these
    Vrec = zeros(pool_options.num_cells,num_timepoints-1);
    Drec_fast = zeros(size(Vrec));
    Drec_slow = zeros(size(Vrec));
    Srec = zeros(size(Vrec));
    
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
                    oflag = false;
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
        
        
        %for recording vars
        Vrec(:,timepoint_counter) = V;
        Drec_fast(:,timepoint_counter) = Dfast;
        Drec_slow(:,timepoint_counter) = Dslow;
        Srec(celltype.excit,timepoint_counter) = Sg_Ex;
        Srec(celltype.inhib,timepoint_counter) = Sg_In;
        
        switch options.fastslow_depression
            case 'off' 
                Dslow(:,idx) =  1; %keep constant at 1.
        end
        
        %check timeout for non-switching or non-dominance 
        if state.count >= state.noswitch_timeout || state.no_dom_counter >= state.no_dominance_timeout
            update_logfile(':::Bistability check failure:::',options.output_log)
            TOF = timepoint_counter*timestep;
            update_logfile(sprintf('---no switch/dominance timeout at t=%.2f(s)',TOF),options.output_log)
            oflag = false;return
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
    
    if ~isempty(durations)
        
        %remove the first artificially induced stay state & subsequent switch state
        trim_Bcheck = find(startsWith(durations(:,end),'leave'), 1, 'first');
        durations = durations(trim_Bcheck+1:end,:);
        sim_results{trialidx,1} = durations;
        switch options.ratelim.check
            case 'on'
            sim_results{trialidx,4} = Rcheck;
        end
        
        
        %taken from find_stay_durations()
        timecourse = size(durations);
        timecourse(2) = timecourse(2) + 1;
        timecourse = cell(timecourse);
        timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
        timecourse(:,3) = cellfun(@(x) mod(x,sum(options.stim_pulse)),timecourse(:,3),'UniformOutput',false);
        %current sample's onset, rounding is needed for subsequent operations
        timecourse(:,4) = cellfun(@(x,y) round(x-y,2),timecourse(:,1),timecourse(:,3),'UniformOutput',false);
        samp_onsets = unique(cat(1,timecourse{:,4})); %like unique won't work properly here without rounding
        timecourse(:,4) = cellfun(@(x) find(x==samp_onsets),timecourse(:,4),'UniformOutput',false); %would also break without rounding
        timecourse(:,end) = durations(:,end);
        
        samp_inds = cell2mat(timecourse(:,4));
        for Sidx = 1:numel(samp_onsets) %unique samples
            curr_rec = samp_inds == Sidx;
            curr_rec = timecourse(curr_rec,:);
            stay_states = startsWith(curr_rec(:,end),'stim');
            if sum(stay_states) > 1 %got heem
                oflag = true;
            end
        end
    end
    oflag = true; %take anything it's fine
    

    %-----unneeded for diagnositcs
    %
    %     %this needs to be fixed, verified like in find_stay_durations()
    %     %check if we need to record a switch index for spiking data
    %     leave_states = find(startsWith(durations(:,end),'leave'));
    %     for LSidx = 1:numel(leave_states)
    %         leave_durr = durations{leave_states(LSidx),2}; %How long leave-state lasted
    %         leave_start = durations{leave_states(LSidx),1} - leave_durr; %when that leave state started
    %         curr_rec = durations(1:leave_states(LSidx),:); %duration record up to that point
    %         prev_stay = find(startsWith(curr_rec(:,end),'stim'), 1, 'last');
    %         %check for no prior stay-state (2 sequential leaves after 1st artificial stay)
    %         if ~isempty(prev_stay)
    %             prev_stay = curr_rec(prev_stay,:); %stay-state prior to leave-transition
    %             last_duration = prev_stay{2}; %how long that stay-state lasted
    %             last_ended = prev_stay{1}; %when it ended
    %             %stay-state followed by a leave-state within Yms later, and lasted at least Y ms
    %             recwin_check = leave_start-last_ended;
    %             recwin_check = recwin_check <= num_postswitch_samples ...
    %                 &&  (recwin_check+leave_durr) >= num_postswitch_samples;
    %             %if the state-state lasted at least Xms, and meets above criteria
    %             if last_duration > num_preswitch_samples && recwin_check
    %                 %record this state in the state-record so we can pull its spiking data out later
    %                 switch_record = vertcat(switch_record,prev_stay);
    %             end
    %         end
    %     end
    %
    %     %-----unneeded for diagnositcs
    %     %check if no switches met recording criteria, terminate if needed
    %     if isempty(switch_record) | numel(switch_record(:,1)) <= 1,return;end
    %     %now get these spike timecourses and save them
    %     num_switches = numel(switch_record(:,1)); %cannot believe this var name is still free
    %     spiking_output = NaN(pool_options.num_cells,num_switch_samples,num_switches);
    %     for tc_idx = 1:num_switches
    %         record_win = switch_record{tc_idx,1};
    %         record_win = record_win-(num_preswitch_samples-1):record_win+num_postswitch_samples;
    %         spiking_output(:,:,tc_idx) = spikes(:,record_win);
    %     end
    
end


if ~oflag %diagnostic conditions unmet
    return
end



update_logfile('---Simulation complete---',options.output_log)
savename = fullfile(options.save_dir,options.sim_name);
save(savename,'sim_results','options')
save(sprintf('%s_D',savename),'Drec_fast','Drec_slow','options','-v7.3')
save(sprintf('%s_V',savename),'Vrec','options','-v7.3')
save(sprintf('%s_S',savename),'Srec','options','-v7.3')
save(sprintf('%s_spikes',savename),'spikes','options','-v7.3')



function state = init_statevar(celltype,options)
timestep = options.timestep;
%---state tracker-------------
state = struct();
state.stay = logical([1 0]);
state.switch = logical([0 1]);
state.undecided = logical([0 0]);
state.last_leave_end = NaN;
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
state.GPU_mdl = options.GPU_mdl;
state.init_check_Lext = options.init_check_Rext * timestep;
state.init_check_stop = options.init_check_tmax / timestep; %minimum time for ready2go check
state.noswitch_timeout = options.noswitch_timeout / timestep;
state.no_dominance_timeout = options.no_dominance_timeout / timestep;
state.no_dom_counter = 0;
state.sample_Estay_offset = options.sample_Estay_offset / timestep; %new: this right here
%same indexing as GPU version for this 
state.pools2compare = [celltype.pool_stay,celltype.pool_switch];
state.pools2compare = state.pools2compare(celltype.excit,:);


