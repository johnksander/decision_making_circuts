clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.

fig_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/Results/NFSv5_combfigs';
if ~isdir(fig_dir)
    mkdir(fig_dir)
end

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

sims2load = {'NFSv5_Estay','NFSv5_Eswitch','NFSv5_Istay','NFSv5_Iswitch'};
simlabels = {'Exit. stay (-)','Exit. switch (+)','Inhib. stay (+)','Inhib. switch (-)'};
num_sims = numel(sims2load);

spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    curr_data = sim_spikerate(sims2load{loadidx},timestep,rtSpikes);
    spike_data{loadidx} = curr_data(:,:,rtSpikes);
end

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%set up for muiltiple figs
fig_fn = 'switching_Xcorr';
figIDs = {'noise','Sg','D','spikes'};
%base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','E-stay * E-switch'};
%I only did noise & spikes, just index these so I dont have to copy & paste etc..
fixvars = [4];
figIDs = figIDs(fixvars);
legloc = legloc(fixvars);
Yax_labs = Yax_labs(fixvars);

%reset this stuff for the plot
num_vars2record = 1;
rtNoise = NaN;rtSpikes = 1; %just so I don't loose track of matrix inds
var_inds = [rtSpikes]; %this is dumb just go with it, don't wana loose track

%only plot -100ms to +100ms
record_duration = size(cat(1,spike_data{:}),2); %get the duration of recorded timecourses
recorded_switchtime =  250e-3/timestep; %actual switchtime
plotting_window = 1+(recorded_switchtime - (100e-3/timestep)):recorded_switchtime + (100e-3/timestep);
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window
%cut down the data matrix to this window
spike_data = cellfun(@(x) x(:,plotting_window),spike_data,'UniformOutput',false);

for figidx = 1:num_vars2record
    
    simRs = cell(num_sims,1);
    %subplot(num_sims/2,num_sims/2,num_sims)
    figure(1)
    orient landscape
    for simidx = 1:num_sims
        
        subplot(num_sims/2,num_sims/2,simidx)
        
        %get the excitable pool neruons, do Xcorr...
        Estay = celltype.excit & celltype.pool_stay;
        Eswitch = celltype.excit & celltype.pool_switch;
        %Istay = celltype.inhib & celltype.pool_stay;
        %Iswitch = celltype.inhib & celltype.pool_switch;
        
        timecourse_data = spike_data{simidx};
        
        Estay = mean(timecourse_data(Estay,:));
        Estay = zeromean(Estay); %zeromean spikerates
        Eswitch = mean(timecourse_data(Eswitch,:));
        Eswitch = zeromean(Eswitch);
        
        [R,laginds] = xcorr(Estay,Eswitch,'unbiased'); %lagging Eswitch
        %         %only plot up to +25ms lag
        %         maxlag = laginds <= 25e-3 / timestep;
        %         laginds = laginds(maxlag);
        %         R = R(maxlag);
        simRs{simidx} = R; %store for another fig
        plot(laginds,R,'Linewidth',2)
        set(gca,'XLim',[min(laginds)-1 max(laginds)+1]);
        %         Xlim_max = numel(R); %this is really annoying, probably remove this if timecourse changes and it's unneeded.
        %         Tzero =  ((Xlim_max - 1)/ 2) + 1;
        %         %onset_switch = 250e-3/timestep; %add the switching time
        %         Xticks = cell(1,5);
        %         Xticks{1} = 0; Xticks{3} = Tzero; Xticks{end} = Xlim_max; %beginning, middle, & end
        %         Xinterval = 150e-3/timestep; %150 ms about zero lag
        %         Xticks{2} = Tzero - Xinterval; Xticks{4} = Tzero + Xinterval;
        %         set(gca,'Xtick',cell2mat(Xticks));
        %         %Xticks = num2cell(get(gca,'Xtick'));
        %         %Xticks{end} = Xlim_max;
        %         Xlabs = cellfun(@(x) sprintf('%+i',round(((x-Tzero)*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        Xticks = num2cell(get(gca,'Xtick'));
        Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        if figidx == rtNoise %just do for noise..
            Yticks = num2cell(get(gca,'Ytick'));
            Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
            set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
        end
        if simidx == 3 || simidx == 4
            xlabel('E-switch lag offset (ms)')
        end
        ylabel(Yax_labs{figidx})
        title(simlabels{simidx})
        
        %y_range = get(gca,'YLim');
        %set(gca,'YLim',[0 max(y_range)]);
        %x_range = get(gca,'XLim');
        %set(gca,'XLim',[0 Xlim_max]);
        %         %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
        set(gca,'Fontsize',fontsz)
        
    end
    
    %align all y axes
    allax = gca;
    allax = allax.Parent.Children;
    linkaxes(allax,'y');
    
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx}]),'-djpeg')
    close all
    figure(2)
    hold on
    for simidx = 1:num_sims
        plot(laginds,simRs{simidx},'Linewidth',2)
    end
    set(gca,'XLim',[min(laginds)-1 max(laginds)+1]);
    Xticks = num2cell(get(gca,'Xtick'));
    %     Xlim_max = numel(simRs{simidx}); %this is really annoying, probably remove this if timecourse changes and it's unneeded.
    %     Tzero =  ((Xlim_max - 1)/ 2) + 1;
    %     %onset_switch = 250e-3/timestep; %add the switching time
    %     Xticks = cell(1,5);
    %     Xticks{1} = 0; Xticks{3} = Tzero; Xticks{end} = Xlim_max; %beginning, middle, & end
    %     Xinterval = 175e-3/timestep; %150 ms about zero lag
    %     Xticks{2} = Tzero - Xinterval; Xticks{4} = Tzero + Xinterval;
    %     set(gca,'Xtick',cell2mat(Xticks));
    %     %Xticks = num2cell(get(gca,'Xtick'));
    %     %Xticks{end} = Xlim_max;
    %     Xlabs = cellfun(@(x) sprintf('%+i',round(((x-Tzero)*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
    Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
    set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    %set(gca,'XLim',[0 Xlim_max]);
    %y_range = get(gca,'YLim');
    %set(gca,'YLim',[0 max(y_range)]);
    set(gca,'Fontsize',fontsz)
    ylabel(Yax_labs{figidx})
    xlabel('E-switch lag offset (ms)')
    legend(simlabels,'location','southeast')
    print(fullfile(fig_dir,[fig_fn '_onefig_' figIDs{figidx}]),'-djpeg')
end







function avg_timecourse = sim_spikerate(sim_name,timestep,rtSpikes)

config_options.modeltype = 'JK';
config_options.sim_name = sim_name;
options = set_options(config_options);

%get results
output_fn = [options.modeltype '_' options.sim_name '.mat'];
sim_output = load(fullfile(options.save_dir,output_fn));
sim_timecourse_output = sim_output.sim_switch_timecourses; %get switching dynamics data

%concatenate results across runs
summed_timecourses = cellfun(@(x) x{1}, sim_timecourse_output, 'UniformOutput', false);
summed_timecourses = cat(4,summed_timecourses{:}); %multivar, cat by 4th dim
summed_timecourses = sum(summed_timecourses,4); %sum by 4th dim, likewise

stateswich_counts = cellfun(@(x) x{2}, sim_timecourse_output, 'UniformOutput', false);
stateswich_counts = sum(cell2mat(stateswich_counts));

avg_timecourse = summed_timecourses ./ stateswich_counts;

% %spikerate for 2ms bins for each cell
num_binsamps = 2e-3/timestep; %num samples in 2ms
cell_spkrates = avg_timecourse(:,:,rtSpikes);
bin_magic = [numel(cell_spkrates(:,1)), num_binsamps, numel(cell_spkrates(1,:))/num_binsamps]; %set up for a magic trick
cell_spkrates = reshape(cell_spkrates,bin_magic);
cell_spkrates = squeeze(sum(cell_spkrates,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_spkrates = repmat(cell_spkrates,[num_binsamps 1 1]);
avg_timecourse(:,:,rtSpikes) = reshape(cell_spkrates,size(avg_timecourse(:,:,rtSpikes))); %put the rabit back in the hat

end

function muC = zeromean(X)
muC = X - mean(X);
end

