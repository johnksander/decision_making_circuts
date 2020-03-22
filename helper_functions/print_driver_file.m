function print_driver_file(options)
%creates a driver file for stim optimization

if numel(options.trial_stimuli) ~= 1,error('not configured correctly for mixed stimuli');end

basedir = fileparts(options.master_driver); %project basedir 

%get options code from master driver 
opt_tags = '%:{3}start:{3}.*%:{3}end:{3}';
master_text = fileread(options.master_driver);

opt_code = regexp(master_text,opt_tags,'match');
opt_code = regexprep(opt_code{:},'\n','\\n'); %insert newlines 
%opt_code = strrep(opt_code,'''',''''''); %I have no idea why this is uneeded
opt_code = strrep(opt_code,'%','%%'); %comments 


outFN = fullfile(options.batchdir,'driver.m');
if exist(outFN,'file') > 0,delete(outFN);end
f = fopen(outFN,'w');

fprintf(f,'clear\nclc\nformat compact\n'); %header 
fprintf(f,'cd %s\n',basedir); %put this in the right directory 
fprintf(f,'jID = str2num([getenv(''SLURM_JOBID''), getenv(''SLURM_ARRAY_TASK_ID'')]);\n'); %jobID
fprintf(f,'idx = %i;\n',options.jobID); %network index from master driver 
fprintf(f,string(opt_code)); %options specification
fprintf(f,'\n\n');
fprintf(f,'options.save_dir = options.batchdir;\n'); %set savedir to batchdir
fprintf(f,'options.sim_name =sprintf(''ES_batch%%i_%%i'',idx,jID);\n'); %reset simname for output files 
fprintf(f,'options.output_log = fullfile(options.save_dir,sprintf(''output_log_%%i.txt'',jID));\n'); %set output log
%add the stimulus intensity 
Rstim = unique(options.trial_stimuli{1});
fprintf(f,'Rstim = %.4f;\n',Rstim); 
fprintf(f,'options.trial_stimuli{1} = [Rstim,Rstim];\n'); 
fprintf(f,'\n\n');

%code for running model 
fprintf(f,'spikeout_model(options);\n'); 

%and output log cleanup
fprintf(f,'delete(options.output_log);\n'); %set output log


fclose(f);
