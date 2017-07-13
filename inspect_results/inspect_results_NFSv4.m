clear
clc
close all
format compact

%result summaries
fontsz = 12;
trial_hists = 'on';
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..
%recorded vars (from noisycurrent_model()):
%lets record, noise, Sg, D, spikes, (I think that's it?)
%just did noise & spikes this time around
num_vars2record = 2;
%rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
%var_inds = [rtNoise,rtSg,rtD,rtSpikes]; %this is dumb just go with it, don't wana loose track
rtNoise = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
var_inds = [rtNoise,rtSpikes]; %this is dumb just go with it, don't wana loose track

orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];


%specify simulation
%---sim setup-----------------
config_options.modeltype = 'JK';
config_options.sim_name = 'NFSv4_Iswitch';
options = set_options(config_options);

%get results
output_fn = [options.modeltype '_' options.sim_name '.mat'];
sim_output = load(fullfile(options.save_dir,output_fn));
options = sim_output.options; %reset to simulation options
options = reset_options_paths(options); %fix paths if needed
sim_timecourse_output = sim_output.sim_switch_timecourses; %get switching dynamics data
sim_output = sim_output.sim_results;
options.trial_hists = trial_hists;
options.conv2secs = timestep;


%Now... make sure the force switch worked. Take only the trials where this is true 
NFS_onset_min = options.NFS_onset_min / timestep; %minimum state duration for forced switch
NFS_stoppush = NFS_onset_min + (options.NFS_stoppush / timestep); %end of noise push

%find stim A & stim B durations
stimA = cellfun(@(x) x{1},sim_output,'UniformOutput', false); %take only stay durations
%get stim A durations
stimA = cellfun(@(x) x(1:2:end),stimA,'UniformOutput', false);
%REMOVE THE FIRST artifically induced switch
stimA = cellfun(@(x) [NaN;x(2:end)],stimA,'UniformOutput', false); %REMOVE THE FIRST artifically induced switch by NaNing it  
%tell me how many are actually valid 
%(this should match number of recorded switches if options.recording window is set correcty)
stimA = cell2mat(stimA);
valid_forcedswitches = stimA >= NFS_onset_min & stimA <= NFS_stoppush;
lnbreak = '-------------------';
fprintf('\r%s\rvalid forced switches = %i\rtotal switches from tarpet state = %i\rovershoots = %i\r%s\r',...
    lnbreak,sum(valid_forcedswitches),numel(valid_forcedswitches),sum(stimA > NFS_stoppush),lnbreak);

constant_switchtime = options.NFS_recordwindow;
constant_switchtime = constant_switchtime / timestep;
constant_switchtime = stimA(valid_forcedswitches) >= min(constant_switchtime) ...
    & stimA(valid_forcedswitches) <= max(constant_switchtime);

fprintf('\r%s\rConstant-time switches = %i\r%s\r',...
    lnbreak,sum(constant_switchtime),lnbreak);
%keeping everything constant.. find switches that occured within a very narrow window
%valid switch averages
%   4.4490e+03 Estay
%   4.4316e+03 Istay
%   4.2999e+03 Iswitch
%   4.3709e+03 Eswitch


% %----- visualization code -----
% %f = histogram(stimA(stimA <  13575) * timestep);
% f = histogram(stimA * timestep);
% f = histogram((stimA(valid_forcedswitches) * timestep) / 1e-3);
% %f.NumBins = 40;
% hold on
% 
% plot(timestep * repmat(NFS_stoppush,numel([min(f.Parent.YLim):1:max(f.Parent.YLim)]),1),...
%     [min(f.Parent.YLim):1:max(f.Parent.YLim)]','linewidth',2,'Color','r')
% plot(timestep * repmat(NFS_onset_min,numel([min(f.Parent.YLim):1:max(f.Parent.YLim)]),1),...
%     [min(f.Parent.YLim):1:max(f.Parent.YLim)]','linewidth',2,'Color','r')
% 
% hold off
% % %----- visualization code -----

%concatenate results across runs
summed_timecourses = cellfun(@(x) x{1}, sim_timecourse_output, 'UniformOutput', false);
summed_timecourses = cat(4,summed_timecourses{:}); %multivar, cat by 4th dim
summed_timecourses = sum(summed_timecourses,4); %sum by 4th dim, likewise

stateswich_counts = cellfun(@(x) x{2}, sim_timecourse_output, 'UniformOutput', false);
stateswich_counts = sum(cell2mat(stateswich_counts));

avg_timecourse = summed_timecourses ./ stateswich_counts;


%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 150;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [2/3 1/3]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);

% %spikerate for 2ms bins for each cell
num_binsamps = 2e-3/timestep; %num samples in 2ms
cell_spkrates = avg_timecourse(:,:,rtSpikes);
bin_magic = [numel(cell_spkrates(:,1)), num_binsamps, numel(cell_spkrates(1,:))/num_binsamps]; %set up for a magic trick
cell_spkrates = reshape(cell_spkrates,bin_magic);
cell_spkrates = squeeze(sum(cell_spkrates,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_spkrates = repmat(cell_spkrates,[num_binsamps 1 1]);
avg_timecourse(:,:,rtSpikes) = reshape(cell_spkrates,size(avg_timecourse(:,:,rtSpikes))); %put the rabit back in the hat 



%set up for muiltiple figs
fig_fn = 'switching_timecourse';
figIDs = {'noise','Sg','D','spikes'};
base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','spike rate (Hz) per 2ms bin'};
%I only did noise & spikes, just index these so I dont have to copy & paste etc.. 
fixvars = [1 4];
figIDs = figIDs(fixvars);
legloc = legloc(fixvars);
Yax_labs = Yax_labs(fixvars);

%only plot -100ms to +100ms 
record_duration = size(avg_timecourse,2);
recorded_switchtime =  250e-3/timestep; %actual switchtime 
plotting_window = 1+(recorded_switchtime - (100e-3/timestep)):recorded_switchtime + (100e-3/timestep);
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window 
%cut down the data matrix to this window 
avg_timecourse = avg_timecourse(:,plotting_window,:);
keyboard
for figidx = 1:num_vars2record
    
    
    %break it down by celltype & aggregate
    Estay = celltype.excit & celltype.pool_stay;
    Eswitch = celltype.excit & celltype.pool_switch;
    Istay = celltype.inhib & celltype.pool_stay;
    Iswitch = celltype.inhib & celltype.pool_switch;
    
    timecourse_data = avg_timecourse(:,:,var_inds(figidx));
    
    Estay = mean(timecourse_data(Estay,:));
    Eswitch = mean(timecourse_data(Eswitch,:));
    Istay = mean(timecourse_data(Istay,:));
    Iswitch = mean(timecourse_data(Iswitch,:));
    %plot the aggregated timecourses
    figure(1)
    hold on
    plot(Estay,'Linewidth',2)
    plot(Eswitch,'Linewidth',2)
    plot(Istay,'Linewidth',2)
    plot(Iswitch,'Linewidth',2)
    hold off
    
    legend({'Excit. stay','Excit. switch','Inhib. stay','Inhib. switch'},'location',legloc{figidx},...
        'Fontsize',fontsz)
    
    Xticks = num2cell(get(gca,'Xtick'));
    Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
    set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    
    if figidx == rtNoise %just do for noise..
        Yticks = num2cell(get(gca,'Ytick'));
        Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
        set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
    end
    xlabel('state switch timecourse (ms)')
    ylabel(Yax_labs{figidx})
    title(base_title)
    y_range = get(gca,'YLim');
    x_range = get(gca,'XLim');
    text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
    set(gca,'Fontsize',fontsz)
    print(fullfile(options.save_dir,[fig_fn '_' figIDs{figidx}]),'-djpeg')
    
    close all
    %now with a smoothed line (only noise)
    %plot the aggregated timecourses
    if figidx == rtNoise
        figure(2)
        hold on
        smoothfac = 10;
        plot(smooth(Estay,smoothfac),'Linewidth',2)
        hold on
        plot(smooth(Eswitch,smoothfac),'Linewidth',2)
        plot(smooth(Istay,smoothfac),'Linewidth',2)
        plot(smooth(Iswitch,smoothfac),'Linewidth',2)
        hold off
        legend({'Excit. stay','Excit. switch','Inhib. stay','Inhib. switch'},'location',legloc{figidx},...
            'Fontsize',fontsz)
        
        Xticks = num2cell(get(gca,'Xtick'));
        Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false);
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        if figidx == rtNoise %just do for noise..
            Yticks = num2cell(get(gca,'Ytick'));
            Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
            set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
        end
        xlabel('state switch timecourse (ms)')
        ylabel(Yax_labs{figidx})
        title(vertcat(base_title,{sprintf('plot smoothed: %i-sample boxcar conv',smoothfac)}))
        y_range = get(gca,'YLim');
        x_range = get(gca,'XLim');
        text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
        set(gca,'Fontsize',fontsz)
        print(fullfile(options.save_dir,[fig_fn '_' figIDs{figidx} '_smooth']),'-djpeg')
        close all
    end
    
    %(only for spikerates)
    %plot the spikerate timecourses scaled down by some factor
    if figidx == rtSpikes
        figure(3)
        scalefac = 2;
        hold on
        plot(Estay ./ scalefac,'Linewidth',2)
        plot(Eswitch ./ scalefac,'Linewidth',2)
        plot(Istay ./ scalefac,'Linewidth',2)
        plot(Iswitch ./ scalefac,'Linewidth',2)
        hold off
        legend({'Excit. stay','Excit. switch','Inhib. stay','Inhib. switch'},'location',legloc{figidx},...
            'Fontsize',fontsz)
        
        Xticks = num2cell(get(gca,'Xtick'));
        Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
        
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        xlabel('state switch timecourse (ms)')
        ylabel(Yax_labs{figidx})
        title(vertcat(base_title,{sprintf('plot scaled down: factor of %i',scalefac)}))
        y_range = get(gca,'YLim');
        x_range = get(gca,'XLim');
        text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
        set(gca,'Fontsize',fontsz)
        print(fullfile(options.save_dir,[fig_fn '_' figIDs{figidx} '_scaled']),'-djpeg')
        close all
        
    end
end


