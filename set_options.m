function options = set_options(config_options)
%basic options template

location = 'bender';

switch location
    case 'lab_desk'
        basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
        options.rand_info = 'shuffle';
    case 'bender'
        basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
        options.rand_info = 'shuffle';
    case 'hpc'
        basedir = '/data/netapp/jksander/rotation/Simulation';
        addpath(fullfile(basedir,'helper_functions')) %annoying...
        options.rand_info = get_rng_seed();
end

rng(options.rand_info)

options.modeltype = config_options.modeltype;
options.sim_name = config_options.sim_name;
options.helper_funcdir = fullfile(basedir,'helper_functions');
results_dir = fullfile(basedir,'Results');
options.save_dir = fullfile(results_dir,options.sim_name);

addpath(options.helper_funcdir)
if ~isdir(options.save_dir),mkdir(options.save_dir);end

%for parameter sweep jobs
if strcmp(options.modeltype,'PS')
    options.jobID = config_options.jobID;
    options.sim_name = sprintf('PS_%s_%i',options.sim_name,options.jobID);
    options.save_dir = fullfile(options.save_dir,options.sim_name);
    if ~isdir(options.save_dir),mkdir(options.save_dir);end
    options.output_log = fullfile(options.save_dir,'output_log.txt');
    
    options.noswitch_timeout = 750; %timeout without a switch (s)
    
    %-----pick connection params-----
    dealers_choice = @(a,b) (a + (b-a).*rand(1));
    
    options.EtoE = .0405; %fixed
    
    options.ItoE = dealers_choice(0.15, 0.65 *2);  %double "fastswitch"
    
    options.EtoI = dealers_choice(0.15, 0.35 *2);  %double "slowswitch"
end

%trial simulation time
if isfield(config_options,'tmax')
    options.tmax = config_options.tmax;
else
    %trial simulation time (s)
    options.tmax = 5000;
end

%whether to force switch from stay state (default false)
if isfield(config_options,'force_back2stay')
    options.force_back2stay = config_options.force_back2stay;
else
    options.force_back2stay = false;
end

%pulse for checking bistability @ sim outset
if isfield(config_options,'init_check_Rext')
    options.init_check_Rext = config_options.init_check_Rext;
else
    options.init_check_Rext = 200; %207.5 hz stimulus to Estay cells
end

%pulse length for checking bistability @ sim outset
if isfield(config_options,'init_check_tmax')
    options.init_check_tmax = config_options.init_check_tmax;
else
    options.init_check_tmax = .9; %pulse must keep steady state for (s)
end

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

update_logfile('--------------------------',options.output_log)


