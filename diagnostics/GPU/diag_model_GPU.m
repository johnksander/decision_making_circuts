function [MF,oflag] = diag_model_GPU(options)
%this model func designed to save whole spiking matrix, for short jobs
MF = mfilename; %for backup purposes
options.GPU_mdl = 'on';
gpurng(options.rand_info)

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
W = gpuArray(reorder_weightmat(W,celltype));
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell basics---------------
Erev = gpuArray(0); %reversal potential, excitatory
Irev = gpuArray(-70e-3); %reversal potential, inhibitory
Gg = gpuArray(10e-9); %max conductance microSiemens
%Tsyn = NaN(pool_options.num_cells,1,'gpuArray'); %gating time constant vector
Tsyn_Ex = gpuArray(50e-3); %excitatory gating time constant, ms
Tsyn_In = gpuArray(10e-3); %inhibitory gating time constant, ms
El = gpuArray(-70e-3); %leak potential mV
Ek = gpuArray(-80e-3); %potassium potential mV
Vreset = gpuArray(-80e-3); %reset potential mV
Rm = gpuArray(100e6); %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = gpuArray(100e-12); %cell capacity picofarads
spike_thresh = gpuArray(20e-3); %spike reset threshold (higher than Vth)
%----noisy input---------------
Tau_ext = NaN(pool_options.num_cells,1,'gpuArray'); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 3.5e-3; % was 2e-3
Tau_ext(celltype.inhib) = 2e-3; % was 5e-3;
initGext = gpuArray(10e-9); %noisy conductance initialization value, nano Siemens
deltaGext = gpuArray(1e-9); %increase noisy conducrance, nano Siemens
Rext = gpuArray(1540); %poisson spike train rate for noise, Hz
%----adaptation conductance----
Vth = gpuArray(-50e-3); %ALEIF spike threshold mV
delta_th = gpuArray(2e-3); %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1,'gpuArray');%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = gpuArray(12.5e-9); %increase adaptation conductance, nano Siemens
%----depression----------------
Pr = gpuArray(.1); %release probability
Td.fast = gpuArray(.3);%synaptic depression time constant, seconds (fast)
Td.slow = gpuArray(10);%slow time constant
fD = gpuArray(.05); %availability factor: (max docked vessicles / max pooled vessicles)
switch options.fastslow_depression
    case 'off'
        fD = 0; Td.slow = Td.fast; 
end
%----timecourse----------------
timestep = gpuArray(options.timestep);
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
    V = zeros(pool_options.num_cells,1,'gpuArray');
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
    stimA = options.trial_stimuli(trialidx,1);
    stimB = options.trial_stimuli(trialidx,2);
    stim_info.stimA_lambda = stimA * timestep; %poisson lambda for stimulus conductance
    stim_info.stimB_lambda = stimB * timestep;
    stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
    if all(~isnan(options.stim_pulse))
        stim_info.delivery = 'pulse';
        stim_info.pulse = options.stim_pulse ./ timestep;
        stim_info.sample_schedule = options.stim_schedule;
    else
        stim_info.delivery = 'constant';
    end
    %---noisy conductance---------
    %Gext = NaN([size(V),2],'gpuArray'); %noisy conductance (do I & E input in 3rd D)
    %Gext(:,1,:) = initGext; %initialize at leak conductance
    ext_inds.I = 1;
    ext_inds.E = 2;
    Gext_Ex = zeros(pool_options.num_cells,1,'gpuArray') + initGext; %initialize at leak conductance
    Gext_In = zeros(pool_options.num_cells,1,'gpuArray') + initGext;
    %---adaptation conductance----
    Gsra = zeros(size(V),'gpuArray');
    %---gating & depression-------
    %Sg = zeros(size(V),'gpuArray'); %synaptic gating
    Sg_Ex = zeros(sum(celltype.excit),1,'gpuArray');
    Sg_In = zeros(sum(celltype.inhib),1,'gpuArray');
    Dfast = ones(size(V),'gpuArray'); %synaptic depression: fast
    Dslow = ones(size(V),'gpuArray'); %synaptic depression: slow
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
    
    %need to preallocate these
    Vrec = zeros(pool_options.num_cells,gather(num_timepoints)-1);
    Drec_fast = zeros(size(Vrec));
    Drec_slow = zeros(size(Vrec));
    Srec = zeros(size(Vrec));
    
    
    while timepoint_counter < num_timepoints
        
        timepoint_counter = timepoint_counter+1;
        state.timeidx = timepoint_counter; %just so I don't have to pass a million things...
        %loop equations
        I = (Erev - V).*(W_Ex*Sg_Ex).*Gg ...
            + (Irev - V).*(W_In*Sg_In).*Gg ...
            + (Gext_Ex.*(Erev-V)) ...
            + (Gext_In.*(Irev-V));
        dVdt = ((El - V + (delta_th.*exp((V - Vth)./delta_th)))./Rm)...
            + (Gsra.*(Ek - V)) + I;
        V = ((dVdt./Cm) .* timestep) + V;
        
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
        %         Dfast(:,idx) = Dfast(:,idx-1) + (((Dslow(:,idx-1) - Dfast(:,idx-1))./Td.fast) .* timestep);%fast syn. depression
        %         Dslow(:,idx) = Dslow(:,idx-1) + timestep .* ( ((1 - Dslow(:,idx-1))./Td.slow) ...
        %             - fD.* ((Dslow(:,idx-1) - Dfast(:,idx-1))./Td.fast)  ); %slow vessicle replacement
        %         Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
        
        Gext_Ex =  Gext_Ex - ((Gext_Ex./Tau_ext) .* timestep); %noisy conductance
        Gext_In =  Gext_In - ((Gext_In./Tau_ext) .* timestep); 
        Gsra = Gsra - ((Gsra./Tsra) .* timestep); %adaptation conductance
        Dfast = Dfast + (((Dslow - Dfast)./Td.fast) .* timestep);%fast syn. depression
        Dslow = Dslow + timestep .* ( ((1 - Dslow)./Td.slow) ...
            - fD.* ((Dslow - Dfast)./Td.fast)  ); %slow vessicle replacement
        Sg_Ex = Sg_Ex - ((Sg_Ex./Tsyn_Ex) .* timestep); %synaptic gating
        Sg_In = Sg_In - ((Sg_In./Tsyn_In) .* timestep);
        
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
            Gext_Ex(celltype.excit) = Gext_Ex(celltype.excit) + (deltaGext.*Pspikes);
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
        Gext_Ex = Gext_Ex + (deltaGext.*ext_spikes(:,ext_inds.E));
        Gext_In = Gext_In + (deltaGext.*ext_spikes(:,ext_inds.I));
        
        
        %for recording vars
        Vrec(:,timepoint_counter) = gather(V);
        Drec_fast(:,timepoint_counter) = gather(Dfast);
        Drec_slow(:,timepoint_counter) = gather(Dslow);
        Srec(celltype.excit,timepoint_counter) = gather(Sg_Ex);
        Srec(celltype.inhib,timepoint_counter) = gather(Sg_In);
        
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
                if timepoint_counter == options.ratelim.stop / timestep
                    [options,Rcheck] = check_rate_limit(spikes,celltype,options);
                    switch Rcheck.status
                        case 'fail'
                            oflag = false;return
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




