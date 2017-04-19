function options = set_options(config_options)
%basic options template 

rng('shuffle')

%basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
%basedir = '/data/netapp/jksander/rotation/Simulation';

options.modeltype = config_options.modeltype;
options.sim_name = config_options.sim_name;
options.helper_funcdir = fullfile(basedir,'helper_functions');
results_dir = fullfile(basedir,'Results');
options.save_dir = fullfile(results_dir,config_options.sim_name);

addpath(options.helper_funcdir)
if ~isdir(options.save_dir)
    mkdir(options.save_dir)
end
