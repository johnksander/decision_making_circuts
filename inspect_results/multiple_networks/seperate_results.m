clear;clc
format compact
close all

%Do you want to seperate simulation results according to some criteria?
%Then this is the place for you. 



%this must be parfored... too many files

basedir = '~/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

num_workers = 24;
c = parcluster('local');
c.NumWorkers = num_workers;
parpool(c,c.NumWorkers,'IdleTimeout',Inf) %,'AttachedFiles',{which('find_stay_durations')})
special_progress_tracker = fullfile(basedir,'inspect_results','multiple_networks','SPT.txt');
if exist(special_progress_tracker) > 0, delete(special_progress_tracker);end %fresh start


sim_name = 'nets_mixstim';
seperate_into = 'netpair-%i'; %which slow/fast network

resdir = fullfile(basedir,'Results',sim_name);
output_fns = dir(fullfile(resdir,'*.mat')); %use this for unrestricted loading
output_fns = cellfun(@(x,y) fullfile(x,y),{output_fns.folder},{output_fns.name},'UniformOutput',false);

num_files = numel(output_fns);
fprintf('-----starting simulation %s (%i files)\n',sim_name,num_files)


gen_options = load(output_fns{1});
gen_options = gen_options.options;
gen_options = rmfield(gen_options,{'stim_targs','trial_stimuli'});
num_net_types = 10;
num_pairs = 5;
pair_inds = num2cell(reshape(1:num_net_types,[],num_pairs)); %just gives a cell array for pair indicies
network_pair_info = cell(num_pairs,1);
for idx = 1:num_pairs
    Pinds = pair_inds(:,idx);
    curr_params = cellfun(@(x) get_network_params(x,gen_options),Pinds,'UniformOutput',false);
    curr_params = cellfun(@(x,y) [y,idx,x.ItoE, x.EtoI],curr_params,Pinds,'UniformOutput',false); 
    curr_params = cat(1,curr_params{:});
    T = array2table(curr_params,'VariableNames',{'net','pair','ItoE','EtoI'});
    network_pair_info{idx} = T;
end
net_info = cat(1,network_pair_info{:});

net_params = net_info{:,{'ItoE','EtoI'}}; 
net_pairs = net_info.pair;

parfor idx = 1:num_files
    
    FN = output_fns{idx};
    curr_file = load(FN);
    
    
    %----seperation criteria
    Fnet = curr_file.options;
    Fnet = [Fnet.ItoE, Fnet.EtoI];
    Fnet = ismember(net_params,Fnet,'rows');
    if sum(Fnet) == 0
        fprintf('WARNING: cannot find match for \n\t%s\n',FN)
        
    else
        %where it goes
        curr_pair = net_pairs(Fnet);      
        outdir = sprintf('%s_%s',resdir,sprintf(seperate_into,curr_pair));
        if ~isdir(outdir),mkdir(outdir);end
        copyfile(FN,outdir)
    end
    
    progress = worker_progress_tracker(special_progress_tracker);
    if mod(progress,floor(num_files * .1)) == 0 %at ten percent
        progress = (progress / num_files) * 100;
        fprintf('%s ---- %.1f percent complete\n',datestr(now,31),progress);
    end
    
end

delete(special_progress_tracker)

delete(gcp('nocreate'))


%for example do::
%sim_name = 'nets_D2t_pref';
%seperate_into = 'stimB-%i'; %B is what percent of A
% %----seperation criteria
% stims = curr_file.options.trial_stimuli;
% B = stims(2) ./stims(1);
% %where it goes
% if round(B*100) == 100
%     outdir = sprintf('%s_%s',resdir,sprintf(seperate_into,round(B*100)));
%     if ~isdir(outdir),mkdir(outdir);end
%     copyfile(FN,outdir)
% end
    
    



