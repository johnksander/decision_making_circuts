clear 
clc
format compact 


profile on

%my model 
%---setup---------------------
jID = 1;
t = 20; %trial simulation time (s) 
options = set_options('modeltype','PS','comp_location','bender',...
    'sim_name','test_profile','jobID',jID,'tmax',t,'netpair_file','D2t');

%------test with stim found for network #1 
options = get_network_params(1,options);
options.EtoE = .0405; %fixed
%---run-----------------------
modelfile = spikeout_model_lite(options);
%---cleanup-------------------
if isempty(dir(fullfile(options.save_dir,'code4*zip')))
    driverfile = mfilename;
    backup_jobcode(options,driverfile,modelfile)
end
%delete(options.output_log) %no need for these right now
logdir = fullfile(options.save_dir,'logs'); %put them seperately
if ~isdir(logdir),mkdir(logdir);end
movefile(options.output_log,logdir)

profile off
profsave(profile('info'),'PROFILE-lite')