function options = set_options(varargin)
%basic options template

%defaults
options.comp_location = 'woodstock';
options.modeltype = ''; %forcing this to break if not specified 
options.sim_name = 'default_name';
options.jobID = 9999; 
options.GPU_mdl = 'off'; %'off' | 'on', set in the actual model code 
options.netpair_file = NaN; %for loading network pair info file 
options.timestep = .02e-3; %.1 milisecond timestep
options.noswitch_timeout = 500; %timeout without a switch (s)
options.no_dominance_timeout = 1; %timeout if neither or both pools active > X seconds 
options.tmax = 1e3; %trial simulation time (s)
options.cut_leave_state = 100e-3; %after Xms in a leave state, half noise E-switch cells 
options.state_def = 'active_states'; %'active_states' | 'include_undecided'; whether simulation aknowledges "undecided states" 
options.state_test_time = 50e-3; %must be X time above threshold to declare a switch 
options.state_test_thresh = .02; %difference in mean Sg between E-cell pools 
options.record_spiking = 'off'; % 'on' saves spiking data, 'off' gives low-memory sim without spiking data 
options.record_preswitch = 250e-3; %transition data recording: 250ms before switch, 150ms after
options.record_postswitch = 150e-3; 
%----pulse stimulus delivery (more realistic licking)
options.stim_pulse = [NaN, NaN]; %default will be none, NaNs specifies constant stim
%format is [on, off] in seconds. options.stim_pulse = [1, 10] gives 1 second pulse w/ 10 second ISI
options.stim_schedule = 'flexible'; % 'fixed' or 'flexible'. Flexible starts stim schedule @ stay-state on  
options.sample_Estay_offset = 40e-3; %(Pulse stim only) init noise offset Estay-Eswitch at the
%sample availablity onset. This is to kick the network into stay & start sampling 
%----checking bistability @ sim outset 
options.init_check_Rext = 400; %pulse strength (Hz to E-stay)
options.init_check_tmax = 1; %pulse must keep steady state for (s)
%----checking spikerates against limits & outset
options.ratelim_check = 'off'; %'off' | 'on'
options.ratelim_E = 50; %hz maximum
options.ratelim_I = 100;
options.ratelim_tmax = 1.5; %time limit (s) for sustained firing over threshold
options.ratelim_mulim = 1.25; % mulim * rate limit gives the average firing rate limit (e.g. 1.25*ratelim.E)
options.ratelim_start = 11; %begin check (s) into sim, check this against init_check_tmax
options.ratelim_stop = 21; %stop check (s) into sim
%----depression
options.fastslow_depression = 'on'; %'on' | 'off' if off, keep Dslow constant at 1 
%if this is off: Dslow should be held constant at 1, tau slow should equal tau fast, and fD = 0. 

%parse inputs 
if mod(numel(varargin),2) ~= 0
    error('arguments must be name-value pairs')
end

num_args = numel(varargin) / 2;
fnames = varargin(1:2:end);
fvals = varargin(2:2:end);
for idx = 1:num_args
    %check for typo first
    if ~isfield(options,fnames{idx})
        error(sprintf('unknown argument: %s',fnames{idx}))
    end
    options.(fnames{idx}) = fvals{idx};
end

nest_fields = {'ratelim'}; %nest_fields = {'LDA','ballistic','TD_l','SVL','Xmemory'};
%now parse the fields you want nested
for idx = 1:numel(nest_fields)
   Fall = fieldnames(options); 
   Fparent = nest_fields{idx};
   F = Fall(startsWith(Fall,Fparent));
   Fnest = cellfun(@(x) strsplit(x,[Fparent '_']),F,'UniformOutput',false);
   Fnest = cellfun(@(x) x{2},Fnest,'UniformOutput',false);
   for j = 1:numel(Fnest)
      options.(Fparent).(Fnest{j}) = options.(F{j});       
   end
   %now remove the holder fields you just nested 
   options = rmfield(options,F);
end


switch options.comp_location 
    case 'woodstock'
        basedir = '/home/acclab/Desktop/ksander/rotation/project';
        options.rand_info = 'shuffle';
    case 'lab_desk'
        basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
        options.rand_info = 'shuffle';
    case 'bender'
        basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
        options.rand_info = 'shuffle';
    case 'hpc64'
        basedir = '/data/netapp/jksander/rotation/Simulation';
        addpath(fullfile(basedir,'helper_functions')) %annoying...
        options.rand_info = get_rng_seed();
    case 'hpc'
        basedir = '/work/jksander/rotation/Simulation';
        addpath(fullfile(basedir,'helper_functions')) %annoying...
        options.rand_info = get_rng_seed();
end

rng(options.rand_info)

options.helper_funcdir = fullfile(basedir,'helper_functions');
results_dir = fullfile(basedir,'Results');
options.save_dir = fullfile(results_dir,options.sim_name);
addpath(options.helper_funcdir)
if ~isdir(options.save_dir),mkdir(options.save_dir);end
options.output_log = fullfile(options.save_dir,sprintf('output_log_%i.txt',options.jobID));


%for parameter sweep jobs
if strcmp(options.modeltype,'PS')
    options.sim_name = sprintf('PS_%s_%i',options.sim_name,options.jobID);    
    %-----pick connection params-----
    dealers_choice = @(a,b) (a + (b-a).*rand(1));
    
    options.EtoE = .0405; %fixed
    switch options.fastslow_depression
        case 'on'
            %---random 
            %range for slow & fast depression sweeep
            %options.ItoE = dealers_choice(0.1, 12.5); 
            %options.EtoI = dealers_choice(0, .75);
            %---grid search
            Ngrid = 100;
            ItoE = linspace(0.1,12.5,Ngrid);
            EtoI = linspace(0,.75,Ngrid);
            [ItoE,EtoI] = meshgrid(ItoE,EtoI);
            ItoE = ItoE(:); EtoI = EtoI(:);
            %HPCC only lets indicies up to 10k!! 
            Gidx = str2num(getenv('SLURM_ARRAY_TASK_ID')); 
            options.ItoE = ItoE(Gidx);
            options.EtoI = EtoI(Gidx);
            
        otherwise
            %range for fast-only depression sweep
            options.ItoE = dealers_choice(0.01, 3.7);
            options.EtoI = dealers_choice(0.01, .35);
    end
    
    %this mode should always be for baseline/no stimulus
    options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'
    Rstim = 0; %rate for stimulus input spikes
    options.trial_stimuli = [Rstim,Rstim];
end

%for stimulus search jobs (e.g. equating network behaviors)
if strcmp(options.modeltype,'equate_stim')
    options.sim_name = sprintf('ES_%s_%i',options.sim_name,options.jobID);
    options.batchdir = fullfile(options.save_dir,sprintf('batch_%i',options.jobID));
    if ~isdir(options.batchdir),mkdir(options.batchdir);end
    
    %-----set network params-----    
    do_config = mod(options.jobID,10);
    do_config(do_config == 0) = 10; 
    options.EtoE = .0405; %fixed
    %pull ItoE, EtoI, Rstim, and stim cell targets for network ID
    options = get_network_params(do_config,options);
end

%for running specific networks (in get_network_params() )
if strcmp(options.modeltype,'NETS') 
    options.sim_name = sprintf('NETS_%s_%i',options.sim_name,options.jobID);
    
    %-----set network params-----    
    do_config = mod(options.jobID,10);
    do_config(do_config == 0) = 10; 
    options.EtoE = .0405; %fixed
    %pull ItoE, EtoI, Rstim, and stim cell targets for network ID
    options = get_network_params(do_config,options);
end

%for diagnostic jobs 
if strcmp(options.modeltype,'diagnostics')
    options.sim_name = sprintf('%s_%i',options.sim_name,options.jobID);
    %diagnostic specific defaults
    if ~ismember('record_spiking',fnames)
        options.record_spiking = 'on'; %default if no argument specified 
    end
end


if sum(strcmp(options.modeltype,{'JK','diagnostics','equate_stim'})) == 0 %don't run this block for these jobs 
    update_logfile('initializing job params...',options.output_log)
    update_logfile(sprintf('tmax = %i',options.tmax),options.output_log)
    update_logfile(sprintf('cut leave after %ims',options.cut_leave_state.*1e3),options.output_log)
    update_logfile('bistability check:',options.output_log)
    update_logfile(sprintf('---pulse = %.1f Hz',options.init_check_Rext),options.output_log)
    update_logfile(sprintf('---test duration = %.1f s',options.init_check_tmax),options.output_log)
    update_logfile(sprintf('---no switch timeout = %i s',options.noswitch_timeout),options.output_log)
    update_logfile('network parameters:',options.output_log)
    if ischar(options.rand_info),message = '---chosen via rng(shuffle)';
    else,message = sprintf('---chosen via rng() seed = %.5f',options.rand_info);end
    update_logfile(message,options.output_log)
    
    update_logfile(sprintf('---EtoE = %.3f',options.EtoE),options.output_log)
    update_logfile(sprintf('---ItoE = %.3f',options.ItoE),options.output_log)
    update_logfile(sprintf('---EtoI = %.3f',options.EtoI),options.output_log)
    if isfield(options,'trial_stimuli')
        update_logfile(sprintf('---trial stimuli = %.1f Hz, %.1f Hz',options.trial_stimuli),options.output_log)
    end
    update_logfile('--------------------------',options.output_log)
    
end

if options.ratelim.start < options.init_check_tmax + 8
    error('Rate limit check is too close to bistability check stimulus. Check slow D tau & relevant parameters.')
    %begin check (s) into sim, check this against init_check_tmax
end

