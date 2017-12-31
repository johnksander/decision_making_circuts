function modelfile = slowswitch_model(options)

modelfile = mfilename; %for backup purposes
options.sim_name = ['JK_' options.sim_name];
output_log = fullfile(options.save_dir,'JK_output_log.txt');
special_progress_tracker = fullfile(options.save_dir,'JK_SPT.txt');
dump_dir = fullfile(options.save_dir,'timecourse_dump'); %only need this for dumping switch timecourses to text files
mkdir(dump_dir)
%make addition stimuli subdirs for dumpdir 
mkdir(fullfile(dump_dir,'stimA'))
mkdir(fullfile(dump_dir,'stimB'))


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
W(W > 0 & EtoE) = (.0675 * .6);
W(W > 0 & ItoE) = (4.15 * .075);
W(W > 0 & EtoI) = .35;
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
tmax = options.tmax; %simulation end (s)
timestep = .25e-3; %.25 milisecond timestep
timevec = 0:timestep:tmax;
num_timepoints = numel(timevec);
Lext = Rext * timestep; %poisson lambda for noisy conductance

%simulation trial loop
%-------------------------------------------------------------------------
update_logfile(':::Starting simulation:::',output_log)
num_trials = numel(options.trial_stimuli(:,1));
sim_results = cell(num_trials,1);

parfor trialidx = 1:num_trials
    
    %preallocate variables
    %-------------------------------------------------------------------------
    %---membrane potential--------
    V = NaN(pool_options.num_cells,2);
    V(:,1) = El; %inital value of membrane potential is leak potential
    %---stimuli info--------------
    stim_info = struct();
    stim_info.targ_cells = celltype.excit & celltype.pool_switch; %Eswitch cells
    stimA = options.trial_stimuli(trialidx,1);
    stimB = options.trial_stimuli(trialidx,2);
    stim_info.stimA_lambda = stimA * timestep; %poisson lambda for stimulus conductance
    stim_info.stimB_lambda = stimB * timestep;
    stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
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
    %disabled for memory usage
    %spikes = zeros(size(V));
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
    state.ready_mintime = .4 / timestep; %minimum time for ready2go check
    state.back2stay_min = .5 / timestep; %wait in a stay state for this long before forcing switch back
    %---switch recording----------
    %250ms before switch, 150ms after
    num_preswitch_samples = 250e-3/timestep;
    num_postswitch_samples = 150e-3/timestep; %also determines how long the push continues after a switch
    num_switch_samples = num_preswitch_samples + num_postswitch_samples;
    %     %lets record, noise, Sg, D, spikes, (I think that's it?)
    %     num_vars2record = 4;
    %     rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
    %only record noise and spikes
    num_vars2record = 2;
    rtInput = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
    stateswich_timecourse = zeros(pool_options.num_cells,num_switch_samples,num_vars2record);
    rolling_timecourse = zeros(size(stateswich_timecourse)); %rolling timecourse matrix
    num_switches_recorded = 0; %for updating the stuff after the switch
    min_NFS_recordwindow = min(options.data_recordwindow / timestep);
    max_NFS_recordwindow = max(options.data_recordwindow / timestep);
    
    experiment_set2go = false; %when experiment is ready to go
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
            %spikes(spiking_cells,idx) = 1;
        end
        
        %update the rolling noise timecourse
        if [timepoint_counter-1] <= num_switch_samples %keep filling it out
            rolling_timecourse(:,timepoint_counter-1,rtInput) = Gext(:,idx);
            %rolling_timecourse(:,timepoint_counter-1,rtSg) = Sg(:,idx);
            %rolling_timecourse(:,timepoint_counter-1,rtD) = D(:,idx);
            rolling_timecourse(:,timepoint_counter-1,rtSpikes) = spiking_cells;
        else
            rolling_timecourse(:,1:end-1,:) = rolling_timecourse(:,2:end,:); %roll it back
            rolling_timecourse(:,end,rtInput) = Gext(:,idx); %add the current noise
            %rolling_timecourse(:,end,rtSg) = Sg(:,idx); %add the current gating
            %rolling_timecourse(:,end,rtD) = D(:,idx); %add the current depression
            rolling_timecourse(:,end,rtSpikes) = spiking_cells; %add the current spikes
        end
        
        [state,durations] = test4switch(Sg(:,idx),state,durations);
        if experiment_set2go
            %find out if we've just switched from A to [the switch before B], Xms ago
            if state.count == num_postswitch_samples & state.now == state.switch
                %if switch Xms ago & we're in the switch state
                prev_state = durations{state.stay}; %find info on the stay state we just switched out of
                last_duration = prev_state{end,1};
                if last_duration >= min_NFS_recordwindow & last_duration <= max_NFS_recordwindow
                    %make sure the switch happened during the recording window
                    prev_stim = prev_state{end,2}; %find the last stimulus type 
                    stateswich_timecourse = stateswich_timecourse + rolling_timecourse; %add the current noise
                    num_switches_recorded = num_switches_recorded + 1; %count it
                    %we're dumping switching timecourses to text files now
                    dumpfn = ['switch_' num2str(trialidx) '_' num2str(num_switches_recorded) '_' num2str(last_duration) '.txt'];
                    dumpfn = fullfile(dump_dir,['stim' prev_stim],dumpfn);
                    dump_data(rolling_timecourse,celltype,rtSpikes,rtInput,dumpfn)
                end
            end
        else
            if state.count >= state.ready_mintime
                durations = cell(1,2); %blank out the recorded durations
                experiment_set2go = true; %we're ready to roll
            end
        end
        
        %input spikes, noise & stimulus
        ext_spikes = poissrnd(Lext,pool_options.num_cells,1); %external spikes to noise conductance
        stim_spikes = timepoint_stimulus(stim_info,state); %get stimulus spikes
        ext_spikes = ext_spikes + stim_spikes; %add 'em both together for one calculation
        ext_spikes = back2stay(ext_spikes,state,celltype); %if in switch state, force back to a stay state
        Gext(:,idx) = Gext(:,idx) + (deltaGext.*ext_spikes); %don't have to index, they get an increase or zero
        
        %lag equation vars for next timepoint
        V = next_timepoint(V);
        Gsra = next_timepoint(Gsra);
        Gext = next_timepoint(Gext);
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


