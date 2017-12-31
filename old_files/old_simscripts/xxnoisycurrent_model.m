function modelfile = noisycurrent_model(options)

modelfile = mfilename; %for backup purposes
options.sim_name = ['JK_' options.sim_name];
output_log = fullfile(options.save_dir,'JK_output_log.txt');
special_progress_tracker = fullfile(options.save_dir,'JK_SPT.txt');
dump_dir = fullfile(options.save_dir,'timecourse_dump'); %only need this for dumping switch timecourses to text files
mkdir(dump_dir)

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
W(W > 0 & EtoE) = (.0675 * .75);   %(.125 / 1.85);
W(W > 0 & ItoE) = (4.15 * .65);    %(20 / 4.8); 
W(W > 0 & EtoI) = (.075 * 2.35);    %(.0375 * 2); 
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
current_val = .1e-9; %applied current in nano amps
noise_sigma = 10e-12; %noise variance in picoamps
current_pulse = options.current_pulse; %switch for adding transient pulses
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
%---adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----timecourse---------------
tmax = options.tmax; %simulation end (s)
timestep = .25e-3; %.25 milisecond timestep
timevec = 0:timestep:tmax;
num_timepoints = numel(timevec);
switch current_pulse
    case 'on'
        update_logfile('WARNING: current pulse functions not configured for low memory',output_log)
end
noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment
update_logfile('WARNING: current_info.noise_sigma set to 0',output_log)
message = sprintf('WARNING: noise sigma at %.1e added outside timepoint_current()',noise_sigma*sqrt(timestep));
update_logfile(message,output_log)
%simulation trial loop
%-------------------------------------------------------------------------
update_logfile(':::Starting simulation:::',output_log)
num_trials = numel(options.trial_currents(:,1));
sim_results = cell(num_trials,1);
sim_switch_timecourses = cell(num_trials,1); %save switching dynamics

parfor trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = NaN(pool_options.num_cells,2);
    V(:,1) = El; %inital value of membrane potential is leak potential
    %---stimuli current ----------
    stimA_curr = options.trial_currents(trialidx,1);
    stimB_curr = options.trial_currents(trialidx,2);
    %---timepoint current info----
    current_info = struct();
    current_info.stimA = stimA_curr;
    current_info.stimB = stimB_curr;
    current_info.basecurr = current_val;
    %BAIT & SWITCH
    %current_info.noise_sigma = noise_sigma; %already adjusted for timestep
    current_info.noise_sigma = 0; %already adjusted for timestep
    %BAIT & SWITCH
    %I'm applying noise outside timepoint_current() so I can see what it is
    
    current_info.initpulse_time = 300; %just index number, should be invariant across units (i.e. PMvsJK)
    current_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
    %---adaptation conductance----
    Gsra = NaN(size(V));
    Gsra(:,1) = 0;
    %---gating & depression-------
    Sg = NaN(size(V)); %synaptic gating
    Sg(:,1) = 0; %initalize at zero??
    D = NaN(size(V)); %synaptic depression
    D(:,1) = 1; %initalize at one
    %---spikes--------------------
    %disabled for memory usage
    %spikes = zeros(size(V));
    %---state tracker-------------
    state.stay = logical([1 0]);
    state.switch = logical([0 1]);
    durations = cell(1,2); %record duration times (stay, switch)
    state.now = state.stay; %pick one state to start with, add pulse to that pool to be sure (make sure this is consistent)
    state.count = 0;
    state.pools2compare = [celltype.pool_stay & celltype.excit,...
        celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
    %---switch recording----------
    %250ms before switch, 100ms after
    num_preswitch_samples = 250e-3/timestep;
    num_postswitch_samples = 100e-3/timestep; %also determines how long the push continues after a switch
    num_switch_samples = num_preswitch_samples + num_postswitch_samples;
    %     %lets record, noise, Sg, D, spikes, (I think that's it?)
    %     num_vars2record = 4;
    %     rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
    %only record noise and spikes
    num_vars2record = 2;
    rtNoise = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
    stateswich_timecourse = zeros(pool_options.num_cells,num_switch_samples,num_vars2record);
    rolling_timecourse = zeros(size(stateswich_timecourse)); %rolling timecourse matrix
    num_switches_recorded = 0; %for updating the stuff after the switch
    %---noise forced switch-------
    current_info.NFS = options.NFS; %cells driving forced switch:  Estay | Eswitch | Istay | Iswitch
    current_info.NFS_onset_min = options.NFS_onset_min / timestep; %minimum state duration for forced switch 
    current_info.NFS_stoppush = current_info.NFS_onset_min + (options.NFS_stoppush / timestep); %end of noise push 
    current_info.NFS_noisepush = options.NFS_noisepush; %noise push
    current_info.NFS_afterswitch_push = num_postswitch_samples; %if switch is successful, keep pushing noise for X samples
    min_NFS_recordwindow = min(options.NFS_recordwindow / timestep);
    max_NFS_recordwindow = max(options.NFS_recordwindow / timestep);
    
    timepoint_counter = 1;
    idx = 2; %keep indexing vars with idx fixed at 2
    
    while timepoint_counter <= num_timepoints
        
        timepoint_counter = timepoint_counter+1;
        
        %get current for this timepoint
        current_info.timeidx = timepoint_counter; %just so I don't have to pass a million things...
        I = timepoint_current(options,current_info,state,durations,celltype);
        curr_noise = noise_sigma * randn(current_info.num_cells,1);
        curr_noise = NFSnoise_adjustment(curr_noise,current_info,state,durations,celltype);
        I = I + curr_noise; %add the noise here
        
        
        %loop equations
        I = I + (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
        I = I + (unique(Erev(celltype.inhib)) - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*unique(Gg(celltype.inhib));
        dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm) + (Gsra(:,idx-1).*(Ek-V(:,idx-1))) + I;
        
        V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
        
        Gsra(:,idx) = Gsra(:,idx-1) - ((Gsra(:,idx-1)./Tsra) .* timestep); %adaptation conductance
        D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
        Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
        
        spiking_cells = V(:,idx) > spike_thresh;
        %changed this if-statement slightly from how it usually is, no functional difference (change needed for rolling timecourse)
        if  sum(spiking_cells) > 0
            Gsra(spiking_cells,idx) = Gsra(spiking_cells,idx) + detlaGsra; %adaptation conductance
            Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
                (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
            D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
            V(spiking_cells,idx) = Vreset;
            %spikes(spiking_cells,idx) = 1;
        end
        
        
        %update the rolling noise timecourse
        if [timepoint_counter-1] <= num_switch_samples %keep filling it out
            rolling_timecourse(:,timepoint_counter-1,rtNoise) = curr_noise;
            %rolling_timecourse(:,timepoint_counter-1,rtSg) = Sg(:,idx);
            %rolling_timecourse(:,timepoint_counter-1,rtD) = D(:,idx);
            rolling_timecourse(:,timepoint_counter-1,rtSpikes) = spiking_cells;
        else
            rolling_timecourse(:,1:end-1,:) = rolling_timecourse(:,2:end,:); %roll it back
            rolling_timecourse(:,end,rtNoise) = curr_noise; %add the current noise
            %rolling_timecourse(:,end,rtSg) = Sg(:,idx); %add the current gating
            %rolling_timecourse(:,end,rtD) = D(:,idx); %add the current depression
            rolling_timecourse(:,end,rtSpikes) = spiking_cells; %add the current spikes
        end
        
        
        if timepoint_counter > current_info.initpulse_time
            [state,durations] = test4switch(Sg(:,idx),state,durations);
            %find out if we've just switched from A to [the switch before B], 100ms ago
            if state.count == num_postswitch_samples & state.now == state.switch
                %if switch 100ms ago & we're in the switch state
                num_stay_states = numel(durations{state.stay});
                if num_stay_states > 2 & mod(num_stay_states,2) == 1
                    %last recorded duration was for stim A & skip over first artifically induced A state
                    %we're in the switch after A now
                    last_duration = durations{state.stay};
                    last_duration = last_duration(end);
                    %if last_duration >= current_info.NFS_onset_min & last_duration <= current_info.NFS_stoppush
                        %make sure the switch happened during the noise push window
                    if last_duration >= min_NFS_recordwindow & last_duration <= max_NFS_recordwindow
                        %make sure the switch happened within a very specfic window, so all timecourses are the same
                        stateswich_timecourse = stateswich_timecourse + rolling_timecourse; %add the current noise
                        num_switches_recorded = num_switches_recorded + 1; %count it
                        %we're dumping switching timecourses to text files now 
                        dumpfn = ['switch_' num2str(trialidx) '_' num2str(num_switches_recorded) '_' num2str(last_duration) '.txt']; 
                        dumpfn = fullfile(dump_dir,dumpfn);
                        dump_data(rolling_timecourse,celltype,rtSpikes,rtNoise,dumpfn)
                    end
                end
            end
        end
        
        %lag equation vars for next timepoint
        V = next_timepoint(V);
        Gsra = next_timepoint(Gsra);
        D = next_timepoint(D);
        Sg = next_timepoint(Sg);
        
    end
    
    sim_results{trialidx} = durations;
    %output summed noise & num samples,use for trial aggregation then averaging
    sim_switch_timecourses{trialidx} = {stateswich_timecourse,num_switches_recorded};
    switch options.parforlog %parfor progress tracking
        case 'on'
            txtappend(special_progress_tracker,'1\n')
            SPT_fid = fopen(special_progress_tracker,'r');
            progress = fscanf(SPT_fid,'%i');
            fclose(SPT_fid);
            if mod(sum(progress),floor(num_trials * .05)) == 0 %5 percent
                progress = (sum(progress) /  num_trials) * 100;
                message = sprintf('Simulation %.1f percent complete',progress);
                update_logfile(message,output_log)
            end
    end
end
update_logfile('---Simulation complete---',output_log)
savename = fullfile(options.save_dir,options.sim_name);
save(savename,'sim_switch_timecourses','sim_results','options')


