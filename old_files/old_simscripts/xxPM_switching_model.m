function xxPM_switching_model(options)

options.sim_name = ['PM_' options.sim_name];
output_log = fullfile(options.save_dir,'PM_output_log.txt');
special_progress_tracker = fullfile(options.save_dir,'PM_SPT.txt');
equal_pools = options.equal_pools; %'on' | 'off'  set switch & stay to have equal properties (b)

%retired 12/28/2016, really high memory usage.

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = NaN; %randomized connection weights, not spatial connections
%--------------------------------------------------------------------------
%build the connectivity matrix
%make celltype logicals
celltype = celltype_logicals(pool_options);
%make these vectors into logical matricies
[EtoE,EtoI,ItoE] = connection_logicals(celltype,pool_options.num_cells);
%make connection scheme based off connection matricies
connection_scheme = EtoE | EtoI | ItoE;  %plan is:  EtoE | EtoI | ItoE;
%make synaptic connection matrix
W = rand(pool_options.num_cells); %use randomized connection weights
W(~connection_scheme) = 0; %filter out connections not allowed by scheme
%modify connection weights
W(EtoE) = .025 * W(EtoE);
W(ItoE) = 6 * W(ItoE);
W(EtoI) = .1 * W(EtoI);
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
current_val = NaN(pool_options.num_cells,1);
current_val(celltype.excit) = 100;
current_val(celltype.inhib) = 80;
noise_sigma = 400; %noise variance in picoamps
current_pulse = options.current_pulse; %switch for adding transient pulses
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential, excitatory
Erev(celltype.inhib) = -62; %reversal potential, inhibitory
Gg = NaN(pool_options.num_cells,1); %max connductance
Gg(celltype.excit) = 1; %just setting these to 1, has no impact.
Gg(celltype.inhib) = 1; %just setting these to 1, has no impact.
Pr = NaN(pool_options.num_cells,1); %release probability
Pr(celltype.excit) = .01; %excitatory release probability
Pr(celltype.inhib) = .01; %inhibitory release probability
W(:,celltype.excit) = W(:,celltype.excit)/unique(Pr(celltype.excit)); %modify connection weights by Pr
Td = 5;%synaptic depression time constant, seconds
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50; %excitatory gating time constant, ms
Tsyn(celltype.inhib) = 5; %inhibitory gating time constant, ms
%----cell basics--------------
El = -62; %leak potential mV
Vreset = NaN(pool_options.num_cells,1); %reset potential mV
Vreset(celltype.excit) = -58; %excitatory
Vreset(celltype.inhib) = -56; %inhibitory
Gl = NaN(pool_options.num_cells,1); %leak conductance
Gl(celltype.excit) = 12; %excitatory
Gl(celltype.inhib) = 10; %inhibitory
Rm = 1./Gl; %resistance
Cm = 200; %cell capacity picofarads
spike_thresh = 20; %spike reset threshold (higher than Vth)
%---adaptation current--------
Vth = -50; %ALEIF spike threshold mV
delta_th = 2; %max voltage threshold, mV  (deltaVth in equation)
a = 2; %adaptation conductance, in nano Siemens
b = NaN(pool_options.num_cells,1); %increase adaptation current by max value, nano amps
b(celltype.excit & celltype.pool_stay) = 3; %excitatory- stay pool
switch equal_pools
    case 'on' %stay and switch pool have same properties
        b(celltype.excit & celltype.pool_switch) = 3; %excitatory- switch pool
    case 'off'%switch pool has less inhibition
        b(celltype.excit & celltype.pool_switch) = 0.5; %excitatory- switch pool
end
b(celltype.inhib) = 0; %inhibitory
Tsra = NaN(pool_options.num_cells,1);%adaptation current time constant, ms
Tsra(celltype.excit) = 300; %excitatory
Tsra(celltype.inhib) = 30; %inhibitory
%---timecourse----------------
tmax = options.tmax; %simulation end, what's the unit here?
timestep = .25; %what's the unit here?
timevec = 0:timestep:tmax;
switch current_pulse
    case 'on'
        %figure out how you want to do this if you wana add these in...
        %pulse_params = options.pulse_params;
        pulse_params = {celltype.excit & celltype.pool_switch,-50,0,200}; %inital pulse for first switch
end
noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment

%simulation trial loop
%-------------------------------------------------------------------------
num_trials = numel(options.trial_currents(:,1));
sim_results = cell(num_trials,1);

parfor trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = NaN(pool_options.num_cells,numel(timevec));
    V(:,1) = -80+40*rand(numel(V(:,1)),1); %inital value of membrane potential is random
    %---applied current vector----
    Iapp = zeros(size(V));
    Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise
    %inhibitory neruons always get the same currrent here!!
    Iapp(celltype.inhib,:) = Iapp(celltype.inhib,:) + unique(current_val(celltype.inhib));
    %---stimuli current ----------
    stimA_curr = options.trial_currents(trialidx,1);
    stimB_curr = options.trial_currents(trialidx,2);
    %---current pulses------------
    switch current_pulse
        case 'on'
            Iapp = addpulse(pulse_params,Iapp,timevec); %add pulses
            %modified current is applied to excitatory cells in both pools for baseline test!
            stimcurr_dif = unique(current_val(celltype.excit))...
                - (unique(current_val(celltype.excit)) * stimA_curr); %current_val vector gets ugly here
            clean_initpulse = {celltype.excit,stimcurr_dif,0,200};
            %modified current is applied to excitatory cells in both pools for baseline test!
            Iapp = addpulse(clean_initpulse,Iapp,timevec); %add pulses
    end
    %---adaptation current--------
    Isra = NaN(size(V));
    Isra(:,1) = 40*rand(numel(Isra(:,1)),1); %initialize at random
    %---gating & depression-------
    Sg = NaN(size(V)); %synaptic gating
    Sg(:,1) = 0;
    D = NaN(size(V)); %synaptic depression
    D(:,1) = 1;
    %---spikes--------------------
    spikes = zeros(size(V));
    %---state tracker-------------
    state.stay = logical([1 0]);
    state.switch = logical([0 1]);
    durations = cell(1,2); %record duration times (stay, switch)
    state.now = state.stay; %pick one state to start with, add pulse to that pool to be sure (make sure this is consistent)
    state.count = 0;
    state.pools2compare = [celltype.pool_stay & celltype.excit,...
        celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
    
    for idx = 2:numel(timevec)
        
        %add noise and preloaded inhibitory current (which is always the same)
        I = Iapp(:,idx-1);
        %add stimulus-specific current
        if mod(numel(durations{state.stay}),2) == 0 %stim A (count is one behind, based off already recorded duration..)
            I(celltype.excit) = I(celltype.excit) + (current_val(celltype.excit) * stimA_curr);
        else %stim B
            I(celltype.excit) = I(celltype.excit) + (current_val(celltype.excit) * stimB_curr);
        end
        
        %loop equations
        I = I + (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
        I = I + (unique(Erev(celltype.inhib)) - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*unique(Gg(celltype.inhib));
        dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm) - Isra(:,idx-1) + I;
        V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
        
        dt_Isra = ((a .* (V(:,idx-1) - El)) - Isra(:,idx-1)) ./ Tsra ; %adaptation current rate of change
        Isra(:,idx) = (dt_Isra .* timestep) + Isra(:,idx-1);
        D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
        Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
        
        if  sum(V(:,idx) > spike_thresh) > 0
            spiking_cells = V(:,idx) > spike_thresh;
            Isra(spiking_cells,idx) = Isra(spiking_cells,idx) + b(spiking_cells); %increase adaptation current by max value
            Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
                (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
            Sg(spiking_cells & celltype.inhib,idx) = 1; %synaptic gating for inhibitory cells
            D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
            D(spiking_cells & celltype.inhib,idx) = 1; %synaptic depression for inhibitory cells
            V(spiking_cells,idx) = Vreset(spiking_cells);
            spikes(spiking_cells,idx) = 1;
        end
        
        [state,durations] = test4switch(Sg(:,idx),state,durations);
    end
    sim_results{trialidx} = durations;
    switch options.parforlog %parfor progress tracking
        case 'on'
            txtappend(special_progress_tracker,'1\n')
            SPT_fid = fopen(special_progress_tracker,'r');
            progress = fscanf(SPT_fid,'%i');
            fclose(SPT_fid);
            if mod(sum(progress),floor(num_trials * .05)) == 0 %5 percent
                progress = (sum(progress) /  num_trials) * 100;
                message = sprintf('Stimulation %.1f percent complete',progress);
                update_logfile(message,output_log)
            end
    end
end
savename = fullfile(options.save_dir,options.sim_name);
save(savename,'sim_results','options')



