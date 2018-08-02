clear
clc
format compact
hold off;close all

%this isn't a "version two", it's just plotting the simulation results all
%together 


%specify simulation
%---sim setup-----------------
Snames = {'network_spiking_P1_1','network_spiking_P2_1','network_spiking_P4_1',...
    'network_spiking_P6_1', 'network_spiking_P8_1', 'network_spiking_P10_1',...
    'network_spiking_P150_1'};
Ftype = 'total_time'; %{'decision_timing','total_samples','total_time'};

figdir = 'network_spiking_pulse_figures';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
addpath(fullfile(basedir,'helper_functions'))
figdir = fullfile(basedir,'Results',figdir,'durations');

%load all the saved datafiles from the other network_spiking_results_durations.m

Nsims = numel(Snames);
sim_data = cellfun(@(x) load(fullfile(figdir,'data',sprintf('%s_%s.mat',x,Ftype))),...
    Snames,'UniformOutput',false);
network_pair_info = cellfun(@(x) x.network_pair_info,sim_data,'UniformOutput',false);
sim_data =  cellfun(@(x) x.result_data,sim_data,'UniformOutput',false);

pulse_durs = cellfun(@(x) strsplit(x,'_'),Snames,'UniformOutput',false);
pulse_durs = cellfun(@(x) str2double(strrep(x{3},'P','')),pulse_durs);


outcome_stat = 'med'; %mu or logmu, med or logmed, var

print_anything = 'yes';

timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb

job_params = cellfun(@(x) x(:,2),sim_data,'UniformOutput',false);
for idx = 1:Nsims
    job_params{idx} = cellfun(@(x)...
        [x.ItoE, x.EtoI,unique(x.trial_stimuli),find(strcmpi(x.stim_targs, stimtarg_vals))],...
        job_params{idx},'UniformOutput',false); %matching "network_pair_info" format
   job_params{idx} = cellfun(@(x) [num2cell(x(1:3)),stimtarg_vals{x(4)}], job_params{idx},'UniformOutput',false);
   job_params{idx} = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Estay','fast')], job_params{idx},'UniformOutput',false);
   job_params{idx} = cellfun(@(x) [x(:,1:3),strrep(x(:,4),'Eswitch','slow')], job_params{idx},'UniformOutput',false);

    
end
%job_params is net_type now 


fig_fn = [Ftype '_%s'];

switch Ftype
    case 'total_samples'
        unit_measure = 'samples';
        %fig_fn = sprintf(fig_fn,'total_samples');
    case 'decision_timing'
        unit_measure = 's - onset';
        %fig_fn = sprintf(fig_fn,'decision_timing');
    case 'total_time'
        unit_measure = 's'; %like "seconds" not samples
        %fig_fn = sprintf(fig_fn,'total_time');
end

switch outcome_stat
    case 'mu'
        fig_fn = Ftype;
        %outcome = mu_duration;
        %Zlabel = 'mean duration (s)';
        Zlabel = sprintf('duration (%s)',unit_measure);
        %figdir = fullfile(figdir,'mean_duration');
    case 'logmu'
        %outcome = logmu_dur;
        %Zlabel = 'mean log duration [ log(s) ]';
        Zlabel = sprintf('log(%s) duration',unit_measure);
        %figdir = fullfile(figdir,'logmean_duration');
        fig_fn = [Ftype,'_log'];
    case 'med'
        %outcome = med_duration;
        Zlabel =  sprintf('med. duration (%s)',unit_measure);
        %figdir = fullfile(figdir,'med_duration');
        fig_fn = [Ftype,'_med'];
    case 'logmed'
        %outcome = med_duration;
        Zlabel =  sprintf('med. duration log(%s)',unit_measure);
        %figdir = fullfile(figdir,'med_duration');
        fig_fn = [Ftype,'_medlog'];
    case 'var'
        Zlabel =  sprintf('duration variance (%s)',unit_measure);
        %figdir = fullfile(figdir,'med_duration');
        fig_fn = [Ftype,'_var'];
end


switch outcome_stat
    case {'med','logmed','var'} %for resampling in parallel
        num_workers = 24;
        c = parcluster('local');
        c.NumWorkers = num_workers;
        parpool(c,c.NumWorkers,'IdleTimeout',Inf)
        BSopt = statset('UseParallel',true);
end

%if ~isdir(figdir),mkdir(figdir);end

matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
BLcol = [103 115 122] ./ 255;
alph = .5;

num_types = unique(cellfun(@numel,network_pair_info));
%plt_idx = 0;
for idx = 1:num_types
    
    curr_net_info = network_pair_info{1}{idx};
    
    %this can select all the ones you want... 
%     curr_net_info = cellfun(@(x) cellfun(@(y) isequal(y,curr_net_info),x),...
%         network_pair_info,'UniformOutput',false);
%     curr_net_info = cellfun(@(x,y) x{y},network_pair_info,curr_net_info,'UniformOutput',false);
%     
    for j = 1:2
        
        %plt_idx = plt_idx + 1;
        h(idx) = subplot(5,1,idx);
        hold on
        
        %get the right color
        if strcmpi(curr_net_info{j,end},'slow')
            lcol = matorange;
        elseif strcmpi(curr_net_info{j,end},'fast')
            lcol = matblue;
        end
        
        %find the right results for network set-up        
        curr_data = cellfun(@(x) cellfun(@(y) isequal(y,curr_net_info(j,:)),x),job_params,'UniformOutput',false);
        curr_data = cellfun(@(x,y) x{y,1}, sim_data,curr_data,'UniformOutput',false);
        %curr_data = cat(1,curr_data{:});
        %curr_data = cell2mat(result_data(curr_data,1));
        
        if ~isempty(curr_data) %skip plot if no data...
            switch outcome_stat
                case {'logmu','logmed'}
                    curr_data = cellfun(@(x) log(x(x ~= 0)),curr_data,'UniformOutput',false);
                    %curr_data = curr_data(curr_data ~= 0);%inf errors
                    %curr_data = log(curr_data);
            end
            
            
            switch outcome_stat
                case {'logmu','mu'}
                    Y = cellfun(@mean,curr_data);
                    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),curr_data);
                case {'logmed','med'}
                    Y = cellfun(@median,curr_data);
                    %SE of the median here, bootstrap estimate
                    SEM = cellfun(@(x) bootstrp(10e3,@median,x,'Options',BSopt),curr_data,'UniformOutput',false);
                    SEM = cellfun(@std,SEM);
                case 'var'
                    Y = cellfun(@var,curr_data);
                    %SE of the variance here, bootstrap estimate
                    SEM = cellfun(@(x) bootstrp(10e3,@var,x,'Options',BSopt),curr_data,'UniformOutput',false);
                    SEM = cellfun(@std,SEM);%should really be "SEV"
            end
            
            errorbar(Y,SEM,'Color',lcol,'LineWidth',3)
            
        end
        
        ylabel({sprintf('network #%i',idx);Zlabel})
        
        set(gca,'XTick',1:Nsims)
        set(gca,'XTickLabel',pulse_durs)
         
        hold off
        leg_info = num2cell(curr_net_info(:,[4,3]),2);
        leg_info = cellfun(@(x) sprintf('%s net: %.0fHz stimulus',x{:}),leg_info,'UniformOutput',false);
        legend(leg_info,'location','best')
        
        if idx == num_types
            xlabel('pulse duration (s)')
        end
        
    end
end
switch outcome_stat
    case {'med','logmed','var'} %for resampling in parallel
        delete(gcp('nocreate'))
end

orient tall
%linkaxes(h,'x')
axis tight

switch print_anything
    case 'yes'
        print(fullfile(figdir,fig_fn),'-djpeg')
end



%summary statistics section, ignore for now

%
% %summary stats
% fprintf('\n------------------------\n')
% fprintf('Summary statistics\n')
% fprintf('%s:\n',Zlabel)
% fprintf('------------------------\n\n')
% %outcome = logmu_dur;
%
% for idx = 1:num_types
%     fprintf('\n------------------------\nNetwork #%i\n',idx)
%     curr_net_info = network_pair_info{idx};
%     Xstim = cell(2,1); %for testing difference between stim distributions
%     for j = 1:2
%         fprintf('---type: %s\n',curr_net_info{j,end})
%         %find the right results for network set-up
%         curr_data = cellfun(@(x) isequal(x,curr_net_info(j,:)),net_type,'UniformOutput',false);
%         curr_data = cat(1,curr_data{:});
%         Xstim{j} = result_data{curr_data,1};
%         fprintf('            stimulus = %.2f\n',outcome(curr_data))
%
%         %take control data
%         BLinfo = [curr_net_info(j,1:2), {0,'baseline'}];
%         curr_data = cellfun(@(x) isequal(x,BLinfo),net_type,'UniformOutput',false);
%         curr_data = cat(1,curr_data{:});
%         fprintf('            baseline = %.2f\n',outcome(curr_data))
%     end
%
%     fprintf('\n---hyp. test: mu stim durrations\n')
%     switch outcome_stat
%         case 'logmu'
%             Xstim = cellfun(@log,Xstim,'UniformOutput',false);
%     end
%     [~,pval] = ttest2(Xstim{1},Xstim{2}); %regular old t-test
%     fprintf('t-test p = %.3f\n',pval)
%     [CI,H] = boot_mean_diffs(Xstim{1},Xstim{2},10000);
%     fprintf('bootstrap test: %s\n',H)
%     fprintf('bootstrap CI: %.2f, %.2f\n',CI)
%
%
% end

