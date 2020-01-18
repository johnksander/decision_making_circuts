clear
clc
format compact
hold off;close all

%This plots results from the stimulus equating jobs. 
basedir = '~/Desktop/work/ACClab/rotation/project/';
curr_dir = fullfile(basedir,'inspect_results/inspect_stimeq');
resdir = fullfile(basedir,'Results/equate_D2t_stims');
RFN = dir(fullfile(resdir,'output_log_*.txt'));
RFN = {RFN.name};

pull_data = 'no';

%net IDs, get the number at the end of output_log_i.txt
NID = cellfun(@(x) strsplit(strrep(x,'.txt',''),'_'),RFN,'UniformOutput',false);
NID = cellfun(@(x) str2double(x{end}),NID);

%need to reformat these for the function below... 
for idx = 1:numel(RFN)
    f = fullfile(resdir,RFN{idx});
    fnew = fullfile(curr_dir,['data-',RFN{idx}]);
    switch pull_data
        case 'yes'
            cmd = sprintf('grep target %s > %s',f,fnew);
            system(cmd);
    end
    plot_stims(fnew)
    T = sprintf('network %i',NID(idx));
    title(T,'FontSize',18,'FontWeight','b')
end




% N = 3;
% for idx = 1:Nsystem(command)
%     
%     plot_stims(sprintf('data%i.txt',idx))
%     
% end



function plot_stims(FN)

%data = readtable('data1.txt','Delimiter','|');

data = readtable(FN);

stim = data.Var4;
T = data.Var7;
T = strrep(T,'s','');
T = str2double(T);

figure;
scatter(stim,T,80,'red','filled')
xlabel('stim Hz'); ylabel('time (s)')
set(gca,'FontSize',18,'FontWeight','b')
end