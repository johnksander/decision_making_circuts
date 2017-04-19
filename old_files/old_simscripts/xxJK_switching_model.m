function xxJK_switching_model(options)

options.sim_name = ['JK_' options.sim_name];
output_log = fullfile(options.save_dir,'JK_output_log.txt');
special_progress_tracker = fullfile(options.save_dir,'JK_SPT.txt');
%retired 12/28/2016, really high memory usage.

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 200;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.5 .5]; %proportion excitable % inhibitory
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
W(W > 0 & EtoE) = .25;
W(W > 0 & ItoE) = 5;
W(W > 0 & EtoI) = .01;
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
%current_val = 0.125e-9; %applied current in nano amps
%noise_sigma = 5e-12; %noise variance in picoamps
current_val = .145e-9; %applied current in nano amps
noise_sigma = 5e-12; %noise variance in picoamps
current_pulse =  options.current_pulse; %switch for adding transient pulses
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
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%---adaptation current--------
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
a = 1.25e-9; %adaptation conductance, in nano Siemens
b = NaN(pool_options.num_cells,1); %increase adaptation current by max value, nano amps
b(celltype.excit) = .25e-9; %excitatory
b(celltype.inhib) = .25e-9; %inhibitory
Tsra = NaN(pool_options.num_cells,1);%adaptation current time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
%----timecourse---------------
tmax = options.tmax; %simulation end (s)
timestep = .25e-3; %.01 milisecond timestep
timevec = 0:timestep:tmax;
switch current_pulse
    case 'on'
        %figure out how you want to do this if you wana add these in...
        %pulse_params = options.pulse_params;
        pulse_params = {celltype.pool_switch,-(current_val/2),0,200e-3}; %inital pulse for first switch
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
    V(:,1) = El; %inital value of membrane potential is leak potential
    %---applied current vector----
    Iapp = zeros(size(V));
    Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise
    %---stimuli current ----------
    stimA_curr = options.trial_currents(trialidx,1);
    stimB_curr = options.trial_currents(trialidx,2);
    %---current pulses------------
    switch current_pulse
        case 'on'
            Iapp = addpulse(pulse_params,Iapp,timevec); %add pulses
            %modified current is applied to both pools in baseline test
            stimcurr_dif = current_val - (current_val * stimA_curr);            
            clean_initpulse = {celltype.pool_switch | celltype.pool_stay,...
                stimcurr_dif,0,200e-3};
             %modified current is applied to both pools in baseline test
             Iapp = addpulse(clean_initpulse,Iapp,timevec); %add pulses
    end
    %---adaptation current--------
    Isra = NaN(size(V));
    Isra(:,1) = 0; %initialize at zero??
    %---gating & depression-------
    Sg = NaN(size(V)); %synaptic gating
    Sg(:,1) = 0; %initalize at zero??
    D = NaN(size(V)); %synaptic depression
    D(:,1) = 1; %initalize at one
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
        
        %add noise
        I = Iapp(:,idx-1);
        %add stimulus-specific current
        if mod(numel(durations{state.stay}),2) == 0 %stim A (count is one behind, based off already recorded duration..)
            I = I + (current_val * stimA_curr);
        else %stim B
            I = I + (current_val * stimB_curr);
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
            D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
            V(spiking_cells,idx) = Vreset;
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


