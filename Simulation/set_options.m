function options = set_options(varargin)
%basic options template

%defaults
options.comp_location = 'woodstock';
options.modeltype = ''; %forcing this to break if not specified 
options.sim_name = 'default_name';
options.jobID = 9999; 
options.timestep = .25e-3; %.25 milisecond timestep
options.noswitch_timeout = 750; %timeout without a switch (s)
options.tmax = 5000; %trial simulation time (s)
options.force_back2stay = false; %whether to force switch from stay state (default false)
options.cut_leave_state = 100e-3; %after Xms in a leave state, cut the noise 
options.state_test_time = 50e-3; %must be X time above threshold to declare a switch 
options.state_test_thresh = .02; %difference in mean Sg between E-cell pools 
options.record_spiking = 'off'; % 'on' saves spiking data, 'off' gives low-memory sim without spiking data 
%----pulse stimulus delivery (more realistic licking) 
options.stim_pulse = [NaN, NaN]; %default will be none
options.stim_schedule = 'flexible'; % 'fixed' or 'flexible' only matters for stim-pulse 
%format is [on, off] in seconds. 
%options.stim_pulse = [1, 10] gives 1 second pulse w/ 10 second ISI
options.sample_Estay_offset = 40e-3; %init noise offset Estay-Eswitch at the
%sample availablity onset. This is to kick the network into stay & start sampling 
%----checking bistability @ sim outset 
options.init_check_Rext = 400; %pulse strength (Hz to E-stay)
options.init_check_tmax = .9; %pulse must keep steady state for (s)
%----checking spikerates against limits & outset
options.ratelim_check = 'off'; %'off' | 'on'
options.ratelim_E = 50; %hz maximum
options.ratelim_I = 100;
options.ratelim_tmax = 1.5; %time limit (s) for sustained firing over threshold
options.ratelim_mulim = 1.25; % mulim * rate limit gives the average firing rate limit (e.g. 1.25*ratelim.E)
options.ratelim_start = 10; %begin check (s) into sim, check this against init_check_tmax
options.ratelim_stop = 20; %stop check (s) into sim
%----depression
options.percent_Dslow = 0; %fraction of slow vesicles (.2 gives 20% slow, 80% fast)
%setting value > 0 enables slower depression timescale. 


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
    if options.percent_Dslow > 0
        %range for slow depression sweeep
        options.ItoE = dealers_choice(0.1, 8);
        options.EtoI = dealers_choice(0.1, 8);
    else
        %range for fast depression sweep
        options.ItoE = dealers_choice(0.1, 4.5);
        options.EtoI = dealers_choice(0.1, 1.5);
    end
    
    %this mode should always be for baseline/no stimulus
    options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'
    Rstim = 0; %rate for stimulus input spikes
    options.trial_stimuli = [Rstim,Rstim];
    
    %from first paramter sweep w/ bad noise
    %options.ItoE = dealers_choice(0.15, 0.65 *2);  %double "fastswitch"
    %options.EtoI = dealers_choice(0.15, 0.35 *2);  %double "slowswitch"
end


if strcmp(options.modeltype,'PS_stim')
    options.sim_name = sprintf('PS_%s_%i',options.sim_name,options.jobID);
    
    %-----set network params-----
    do_config = mod(options.jobID,10);
    options.EtoE = .0405; %fixed
    %pull ItoE, EtoI, Rstim, and stim cell targets for network ID 
    options = get_network_params(do_config,options); 
end


if strcmp(options.modeltype,'diagnostics')
    options.sim_name = sprintf('%s_%i',options.sim_name,options.jobID);
    %diagnostic specific defaults
    if ~ismember('record_spiking',fnames)
        options.record_spiking = 'on'; %default if no argument specified 
    end
end


if sum(strcmp(options.modeltype,{'JK','diagnostics'})) == 0 %don't run this block for outdated jobs 
    update_logfile('initializing job params...',options.output_log)
    update_logfile(sprintf('tmax = %i',options.tmax),options.output_log)
    update_logfile(sprintf('force back2stay = %s',string(options.force_back2stay)),options.output_log)
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

