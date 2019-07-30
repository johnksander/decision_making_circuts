clear
clc
format compact



Snames = {'nets_D2t_pref'}; %Snames = {'nets_fastD','nets_slowD'};

basedir = '/home/acclab/Desktop/ksander/rotation/project';

for idx = 1:numel(Snames)
    fix_this(basedir,Snames{idx})
end


function fix_this(basedir,sim_name)
resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,['*',sim_name,'*.mat'])); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);
BL_fns = dir(fullfile([resdir '_baseline'],['*',sim_name,'*.mat']));
BL_fns = cellfun(@(x,y) fullfile(x,y),{BL_fns.folder},{BL_fns.name},'UniformOutput',false);
output_fns = cat(2,BL_fns,output_fns);

fprintf('-----starting simulation results: %s\n',sim_name)

num_files = numel(output_fns);
num_deleted = 0;
for idx = 1:num_files
    if mod(idx,500) == 0,fprintf('working on file #%i/%i...\n',idx,num_files);end
    curr_file = load(output_fns{idx});
    %get state durations
    try
        state_durations = curr_file.sim_results;
    catch
        
        delete(output_fns{idx})
        num_deleted = num_deleted + 1;
    end
    
end

fprintf('\n\n-----total deleted: %i out of %i total files\n\n',num_deleted,num_files)

end
