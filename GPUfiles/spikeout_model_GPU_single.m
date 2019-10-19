function modelfile = spikeout_model_GPU_single(options)
%this model func designed to save whole spiking matrix, for short jobs
modelfile = mfilename; %for backup purposes
options.GPU_mdl = 'on';
gpurng(options.rand_info)

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
W = gpuArray(reorder_weightmat(W,celltype));
%celltype = structfun(@(x) gpuArray(x),celltype,'UniformOutput',false);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell basics---------------
Zvec = zeros(pool_options.num_cells,1,'single','gpuArray'); %for GPU calcs
Zvec_Ex = zeros(sum(celltype.excit),1,'single','gpuArray');
Zvec_In = zeros(sum(celltype.inhib),1,'single','gpuArray'); 
Erev = 0 + Zvec; %reversal potential, excitatory
Irev = -70e-3 + Zvec; %reversal potential, inhibitory
Gg = 10e-9 + Zvec; %max conductance microSiemens
%Tsyn = NaN(pool_options.num_cells,1,'single','gpuArray'); %gating time constant vector
Tsyn_Ex = 50e-3 + Zvec_Ex; %excitatory gating time constant, ms
Tsyn_In = 10e-3 + Zvec_In; %inhibitory gating time constant, ms
El = -70e-3 + Zvec; %leak potential mV
Ek = -80e-3 + Zvec; %potassium potential mV
Vreset = gpuArray(-80e-3); %reset potential mV
Rm = 100e6 + Zvec; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12 + Zvec; %cell capacity picofarads
spike_thresh = gpuArray(20e-3); %spike reset threshold (higher than Vth)
%----noisy input---------------
Tau_ext = NaN(pool_options.num_cells,1,'single','gpuArray'); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 3.5e-3; % was 2e-3
Tau_ext(celltype.inhib) = 2e-3; % was 5e-3;
initGext = gpuArray(10e-9); %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9 + Zvec; %increase noisy conducrance, nano Siemens
Rext = gpuArray(1540); %poisson spike train rate for noise, Hz
%----adaptation conductance----
Vth = -50e-3 + Zvec; %ALEIF spike threshold mV
delta_th = 2e-3 + Zvec; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1,'single','gpuArray');%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = gpuArray(12.5e-9); %increase adaptation conductance, nano Siemens
%----depression----------------
Pr = gpuArray(.1); %release probability
Td_fast = .3 + Zvec;%synaptic depression time constant, seconds (fast)
Td_slow = 10 + Zvec;%slow time constant
fD = .05 + Zvec; %availability factor: (max docked vessicles / max pooled vessicles)
switch options.fastslow_depression
    case 'off'
        fD = 0; Td_slow = Td_fast;
end
%----timecourse----------------
timestep = gpuArray(options.timestep);
dtvec = timestep + Zvec; %vectors for GPU calcs
dt_Ex = timestep + Zvec_Ex;
dt_In = timestep + Zvec_In;
num_timepoints = gather(round(options.tmax / timestep) + 1); %same as numel(0:dt:Tmax)
Lext = Rext * timestep; %poisson lambda for noisy conductance
check_GPU_lambda(Lext); %check lambda 
Lext = repmat(Lext,pool_options.num_cells,2);

%simulation trial loop
%-------------------------------------------------------------------------
update_logfile(':::Starting simulation:::',options.output_log)
num_trials = numel(options.trial_stimuli(:,1));
sim_results = cell(num_trials,4);

for trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = zeros(pool_options.num_cells,1,'single','gpuArray');
    V = V + El;%inital value of membrane potential is leak potential
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
    stim_info.targ_cells = gpuArray(stim_info.targ_cells);
    stimA = options.trial_stimuli(trialidx,1) * timestep; %poisson lambda for stimulus conductance
    stimB = options.trial_stimuli(trialidx,2) * timestep;
    check_GPU_lambda(stimA); check_GPU_lambda(stimB); %check lambda 
    stim_info.stimA_lambda = repmat(stimA,sum(stim_info.targ_cells),1);
    stim_info.stimB_lambda = repmat(stimB,sum(stim_info.targ_cells),1);
    stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
    if all(~isnan(options.stim_pulse))
        stim_info.delivery = 'pulse';
        stim_info.pulse = options.stim_pulse ./ timestep;
        stim_info.sample_schedule = options.stim_schedule;
    else
        stim_info.delivery = 'constant';
    end
    %---noisy conductance---------
    %Gext = NaN([size(V),2],'single','gpuArray'); %noisy conductance (do I & E input in 3rd D)
    %Gext(:,1,:) = initGext; %initialize at leak conductance
    ext_inds.I = 1;
    ext_inds.E = 2;
    Gext_Ex = initGext + Zvec; %initialize at leak conductance
    Gext_In = initGext + Zvec;
    %---adaptation conductance----
    Gsra = zeros(size(V),'single','gpuArray');
    %---gating & depression-------
    %Sg = zeros(size(V),'single','gpuArray'); %synaptic gating
    Sg_Ex = zeros(sum(celltype.excit),1,'single','gpuArray');
    Sg_In = zeros(sum(celltype.inhib),1,'single','gpuArray');
    Dfast = ones(size(V),'single','gpuArray'); %synaptic depression: fast
    Dslow = ones(size(V),'single','gpuArray'); %synaptic depression: slow
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
    %idx = 2; %keep indexing vars with idx fixed at 2
    
    %equation vars (move these when done)
    W_Ex = W(:,celltype.excit);
    W_In = W(:,celltype.inhib);
    while timepoint_counter < num_timepoints
        
        timepoint_counter = timepoint_counter+1;
        state.timeidx = timepoint_counter; %just so I don't have to pass a million things...
        %loop equations
        
        WS_Ex = W_Ex*Sg_Ex;
        WS_In = W_In*Sg_In;
        I = arrayfun(@IFUN,V,Erev,Irev,WS_Ex,WS_In,Gg,Gext_Ex,Gext_In);
        V = arrayfun(@VtFUN,El,V,delta_th,Vth,Rm,Gsra,Ek,I,Cm,dtvec);
        
        Gext_Ex = arrayfun(@XtFUN,Gext_Ex,Tau_ext,dtvec);%noisy conductance
        Gext_In = arrayfun(@XtFUN,Gext_In,Tau_ext,dtvec);
        Gsra = arrayfun(@XtFUN,Gsra,Tsra,dtvec);%adaptation conductance
        Sg_Ex = arrayfun(@XtFUN,Sg_Ex,Tsyn_Ex,dt_Ex);%synaptic gating
        Sg_In = arrayfun(@XtFUN,Sg_In,Tsyn_In,dt_In);
        [Dfast,Dslow] = arrayfun(@DtFUN,Dfast,Dslow,Td_fast,Td_slow,fD,dtvec);%fast & slow syn depression
        
        %         I = (Erev - V).*(W_Ex*Sg_Ex).*Gg ...
        %             + (Irev - V).*(W_In*Sg_In).*Gg ...
        %             + (Gext_Ex.*(Erev-V)) ...
        %             + (Gext_In.*(Irev-V));
        %         dVdt = ((El - V + (delta_th.*exp((V - Vth)./delta_th)))./Rm)...
        %             + (Gsra.*(Ek - V)) + I;
        %         V = ((dVdt./Cm) .* timestep) + V;
        
        %         I = (Erev - V(:,idx-1)).*(W_Ex*Sg(celltype.excit,idx-1)).*Gg;
        %         I = I + (Irev - V(:,idx-1)).*(W_In*Sg(celltype.inhib,idx-1)).*Gg;
        %         I = I + (Gext(:,idx-1,ext_inds.E).*(Erev-V(:,idx-1)));
        %         I = I + (Gext(:,idx-1,ext_inds.I).*(Irev-V(:,idx-1)));
        %         dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm)...
        %             + (Gsra(:,idx-1).*(Ek-V(:,idx-1))) + I;
        %
        %         V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
        
        
        %         tGext = squeeze(Gext(:,idx-1,:)); %flattened t-1 Gext
        %         Gext(:,idx,:) = tGext - ((tGext./Tau_ext) .* timestep); %noisy conductance
        %         Gsra(:,idx) = Gsra(:,idx-1) - ((Gsra(:,idx-1)./Tsra) .* timestep); %adaptation conductance
        %         Dfast(:,idx) = Dfast(:,idx-1) + (((Dslow(:,idx-1) - Dfast(:,idx-1))./Td_fast) .* timestep);%fast syn. depression
        %         Dslow(:,idx) = Dslow(:,idx-1) + timestep .* ( ((1 - Dslow(:,idx-1))./Td_slow) ...
        %             - fD.* ((Dslow(:,idx-1) - Dfast(:,idx-1))./Td_fast)  ); %slow vessicle replacement
        %         Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
        
        %         Gext_Ex =  Gext_Ex - ((Gext_Ex./Tau_ext) .* timestep); %noisy conductance
        %         Gext_In =  Gext_In - ((Gext_In./Tau_ext) .* timestep);
        %         Gsra = Gsra - ((Gsra./Tsra) .* timestep); %adaptation conductance
        %         Dfast = Dfast + (((Dslow - Dfast)./Td_fast) .* timestep);%fast syn. depression
        %         Dslow = Dslow + timestep .* ( ((1 - Dslow)./Td_slow) ...
        %             - fD.* ((Dslow - Dfast)./Td_fast)  ); %slow vessicle replacement
        %         Sg_Ex = Sg_Ex - ((Sg_Ex./Tsyn_Ex) .* timestep); %synaptic gating
        %         Sg_In = Sg_In - ((Sg_In./Tsyn_In) .* timestep);
        
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
            [state,durations] = test4switch(gather(Sg_Ex),state,durations);
            [state,avail_noise] = check_noise_avail(stim_info,state);
        else
            state.count = state.count + 1; %if you don't run test4switch(), must update this counter outside
        end
        
        %run the bistability check
        if ~experiment_set2go %during bistability check, check_bistability() handles pulse input spikes
            [BScheck,Pspikes,state] = check_bistability(Sg_Ex,state);
            %add pulse spikes (same as below)
            Gext_Ex(celltype.excit) = Gext_Ex(celltype.excit) + unique(deltaGext).*Pspikes;
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
        ext_spikes = arrayfun(@POISFUN,Lext);
        
        %adjust noise input if needed
        if avail_noise.Estay ~= 1
            ext_spikes(celltype.excit & celltype.pool_stay,:) = ...
               arrayfun(@POISFUN,avail_noise.Estay.*Lext);
        end
        if avail_noise.Eswitch ~= 1
            ext_spikes(celltype.excit & celltype.pool_switch,:) = ...
                arrayfun(@POISFUN,avail_noise.Eswitch.*Lext);
        end
        
        %---spiking input from stimulus
        if experiment_set2go
            stim_spikes = timepoint_stimulus(stim_info,state); %get stimulus spikes
            %always exitatory, add 'em both together for one calculation
            ext_spikes(:,ext_inds.E) = ext_spikes(:,ext_inds.E) + stim_spikes;
        end
        
        %update Gexternal. Don't have to index, they get an increase or zero
        [Gext_Ex,Gext_In] = arrayfun(@GSPIKEFUN,Gext_Ex,Gext_In,deltaGext,...
            ext_spikes(:,ext_inds.E),ext_spikes(:,ext_inds.I));
        %         Gext_Ex = Gext_Ex + (deltaGext.*ext_spikes(:,ext_inds.E));
        %         Gext_In = Gext_In + (deltaGext.*ext_spikes(:,ext_inds.I));
        
%         %lag equation vars for next timepoint
%         V = next_timepoint(V);
%         Gsra = next_timepoint(Gsra);
%         Gext = next_timepoint(Gext);
%         Dfast = next_timepoint(Dfast);
%         Dslow = next_timepoint(Dslow);
%         Sg = next_timepoint(Sg);
        
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
    if isempty(durations),return;end
    trim_Bcheck = find(startsWith(durations(:,end),'leave'), 1, 'first');
    durations = durations(trim_Bcheck+1:end,:); if isempty(durations),return;end
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
            num_preswitch_samples = options.record_preswitch/timestep; %window samples
            num_postswitch_samples = options.record_postswitch/timestep;
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

function stim_spikes = timepoint_stimulus(stim_info,state)
%inputs: stimuli info & state parameter structures

switch state.GPU_mdl
    case 'off'
        stim_spikes = zeros(stim_info.num_cells,1);
    case 'on'
        stim_spikes = zeros(stim_info.num_cells,1,'single','gpuArray');
end

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
    %add stimulus-specific spiketrain
    if strcmp(state.stim_labels{state.current_stimulus},'stim_A')
        stim_spikes(stim_info.targ_cells) = arrayfun(@POISFUN,stim_info.stimA_lambda);
    elseif strcmp(state.stim_labels{state.current_stimulus},'stim_B')
        stim_spikes(stim_info.targ_cells) = arrayfun(@POISFUN,stim_info.stimB_lambda);
    end
end

function [BScheck,Pspikes,state] = check_bistability(Sg,state)

BScheck.status = 'running'; %initalize
target_cells = state.pools2compare(:,state.stay);
Pspikes = zeros(size(target_cells),'single','gpuArray');
%see if we're still in the bistability check window
check_window = state.count <= state.init_check_stop;
%add another one for "check over", just in case isnan(state.count) or something...
check_over = state.count >= state.init_check_stop;
%give pulse up until 2/3 into check window, then stop
pulse_window = state.count <= floor(state.init_check_stop *8); % *.66 would give 2/3
%if the pulse window is 1/3 over, start checking state dominance
check_dom = state.count >= floor(state.init_check_stop *.33) && check_window;


if pulse_window
    %give pulse spikes
    Pspikes(target_cells) = arrayfun(@POISFUN,state.init_check_Lext); %external spikes to noise conductance
end

if check_dom
    %test if there's an active state, check if it's the targeted cells (E-stay)
    Smu = [Sg(state.pools2compare(:,1)), Sg(state.pools2compare(:,2))];
    Smu = mean(Smu);
    active_state = Smu - fliplr(Smu); %A-B, B-A
    state.now = active_state > state.test_thresh; %will return undecided, if neither active
    
    targ_active = all(state.now == state.stay); %if it's the state we intended
    if ~targ_active
        state.thresh_clock = state.thresh_clock + 1; 
    elseif targ_active
        state.thresh_clock = 0; %reset the clock
    end
end

%network shouldn't leave target state during the check period
no_switches = state.thresh_clock < state.test_time;


%---outcomes---

if check_dom && ~no_switches
    %active network state isn't the targeted cells---- fail
    if ~strcmp(BScheck.status,'fail')
        %make sure hasn't already failed, give correct message
        BScheck.status = 'fail';
        muA = mean(Sg(state.pools2compare(:,state.stay)));
        muB = mean(Sg(state.pools2compare(:,state.switch)));
        BScheck.Fcode = sprintf('---non-dominance (muStay=%.2f,muLeave=%.2f)',muA,muB);
    end
end

%this was an outcome under the old switch-testing scheme 
% if ~no_switches && check_window
%     %active network state isn't the targeted cells during the pulse window---- fail
%     BScheck.status = 'fail';
%     BScheck.Fcode = '---N > 1 switches';
% end

no_failure = ~strcmp(BScheck.status,'fail');

if check_over && no_switches && no_failure
    %should be good, don't forget to renew your sticker next year
    BScheck.status = 'pass';
end

function check_GPU_lambda(lambda)
nope = lambda >= 15 | lambda < 0 | isinf(lambda); %GPU considerations
if any(nope),error('cannot handle lambda');end

function I = IFUN(V,Erev,Irev,WS_Ex,WS_In,Gg,Gext_Ex,Gext_In)
E_diff = Erev-V;
I_diff = Irev-V;
I = E_diff.*WS_Ex.*Gg + I_diff.*WS_In.*Gg + Gext_Ex.*E_diff + Gext_In.*I_diff;

function V = VtFUN(El,V,delta_th,Vth,Rm,Gsra,Ek,I,Cm,dt)

dVdt = ((El - V + (delta_th.*exp((V - Vth)./delta_th)))./Rm)...
    + (Gsra.*(Ek - V)) + I;

V = ((dVdt./Cm) .* dt) + V;

function X = XtFUN(X,tau,dt)
%generic function for updates in the form :
% Gnew = G - (G / tau) * timestep   (e.g. conductance)

X = X - (X./tau) .* dt;

function [Dt_fast,Dt_slow] = DtFUN(Dfast,Dslow,tauD_fast,tauD_slow,fD,dt)

Ddiff = Dslow - Dfast;

Dt_fast = Dfast + ((Ddiff./tauD_fast) .* dt);%fast syn. depression

Dt_slow = Dslow + dt .* ( ((1 - Dslow)./tauD_slow) ...
    - fD.* (Ddiff./tauD_fast)  ); %slow vessicle replacement

function [Gext_Ex,Gext_In] = GSPIKEFUN(Gext_Ex,Gext_In,deltaGext,spikes_Ex,spikes_In)

Gext_Ex = Gext_Ex + deltaGext.*spikes_Ex;
Gext_In = Gext_In + deltaGext.*spikes_In;

function t = POISFUN(lambda)
t = 0;
p = 0 - log(rand());
while p < lambda
    t = t + 1; 
    p = p - log(rand());
end
