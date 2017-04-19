function backup_jobcode(options,driverfile,modelfile)
%creates a zip file of the driver, model, and helper function dir. 
%stores in results folder 

modelfile = which(modelfile);
driverfile = which(driverfile);

zipname = ['code4' options.sim_name];
zipname = fullfile(options.save_dir,zipname);

zip(zipname,{modelfile,driverfile,options.helper_funcdir})




