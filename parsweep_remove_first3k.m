clear
clc
format compact


%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_baseline';
basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
resdir = fullfile(basedir,'Results',sim_name);
sepdir = fullfile(basedir,'Results','parsweep_first3k');
if ~isdir(sepdir),mkdir(sepdir);end

output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat']));
output_fns = {output_fns.name};
job_num = cellfun(@(x) strsplit(x,'PS_parsweep_baseline_'),output_fns,'UniformOutput',false);
job_num = cellfun(@(x) strsplit(x{end},'.mat'),job_num,'UniformOutput',false);
job_num = cellfun(@(x) str2double(x{1}),job_num,'UniformOutput',false);
job_num = cell2mat(job_num);

files2move = job_num < 3001;
files2move = output_fns(files2move);
for idx = 1:numel(files2move)
    
    
    status = movefile(fullfile(resdir,files2move{idx}),fullfile(sepdir,files2move{idx}));
    if status == 0
        fprintf('failure: %s\r',files2move{idx})
    end    
end


