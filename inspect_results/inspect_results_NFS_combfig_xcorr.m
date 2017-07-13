clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.

fig_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/Results/NFSv3_combfigs';
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

sims2load = {'NFSv3_Estay','NFSv3_Eswitch','NFSv3_Iswitch','NFSv3_Istay'};
simlabels = {'Exit. stay (-)','Exit. switch (+)','Inhib. switch (-)','Inhib. stay (+)'};
num_sims = numel(sims2load);

spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    spike_data{loadidx} = sim_spiketrain(sims2load{loadidx},rtSpikes);
    %spike_data{loadidx} = curr_data(:,:,rtSpikes);
end

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 150;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [2/3 1/3]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%set up for muiltiple figs
fig_fn = 'switching_Xcorr';
figIDs = {'noise','Sg','D','spikes'};
%base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','spike count'};
%I only did noise & spikes, just index these so I dont have to copy & paste etc..
fixvars = [4];
figIDs = figIDs(fixvars);
legloc = legloc(fixvars);
Yax_labs = Yax_labs(fixvars);

%reset this stuff for the plot
num_vars2record = 1;
rtNoise = NaN;rtSpikes = 1; %just so I don't loose track of matrix inds
var_inds = [rtSpikes]; %this is dumb just go with it, don't wana loose track


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
        
        Estay = sum(timecourse_data(Estay,:));
        Eswitch = sum(timecourse_data(Eswitch,:));
        R = xcorr(Estay,Eswitch,'coeff');
        %R = xcorr(Estay,Eswitch,'unbiased');
        simRs{simidx} = R; %store for another fig
            
        plot(R,'Linewidth',2)
        
        Xlim_max = numel(R); %this is really annoying, probably remove this if timecourse changes and it's unneeded.
        Tzero =  ((Xlim_max - 1)/ 2) + 1;
        %onset_switch = 250e-3/timestep; %add the switching time
        Xticks = cell(1,5);
        Xticks{1} = 0; Xticks{3} = Tzero; Xticks{end} = Xlim_max; %beginning, middle, & end
        Xinterval = 150e-3/timestep; %150 ms about zero lag
        Xticks{2} = Tzero - Xinterval; Xticks{4} = Tzero + Xinterval;
        set(gca,'Xtick',cell2mat(Xticks));
        %Xticks = num2cell(get(gca,'Xtick'));
        %Xticks{end} = Xlim_max;
        Xlabs = cellfun(@(x) sprintf('%+i',round(((x-Tzero)*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        if figidx == rtNoise %just do for noise..
            Yticks = num2cell(get(gca,'Ytick'));
            Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
            set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
        end
        if simidx == 3 || simidx == 4
            xlabel('lag offset (ms)')
        end
        ylabel(Yax_labs{figidx})
        title(simlabels{simidx})
        
        y_range = get(gca,'YLim');
        set(gca,'YLim',[0 max(y_range)]);
        %x_range = get(gca,'XLim');
        set(gca,'XLim',[0 Xlim_max]);
        %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
        set(gca,'Fontsize',fontsz)

    end
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx}]),'-djpeg')
    close all
    figure(2)
    hold on
    for simidx = 1:num_sims
        plot(simRs{simidx},'Linewidth',2)
    end
    Xticks = num2cell(get(gca,'Xtick'));
    
    Xlim_max = numel(simRs{simidx}); %this is really annoying, probably remove this if timecourse changes and it's unneeded.
    Tzero =  ((Xlim_max - 1)/ 2) + 1;
    %onset_switch = 250e-3/timestep; %add the switching time
    Xticks = cell(1,5);
    Xticks{1} = 0; Xticks{3} = Tzero; Xticks{end} = Xlim_max; %beginning, middle, & end
    Xinterval = 175e-3/timestep; %150 ms about zero lag
    Xticks{2} = Tzero - Xinterval; Xticks{4} = Tzero + Xinterval;
    set(gca,'Xtick',cell2mat(Xticks));
    %Xticks = num2cell(get(gca,'Xtick'));
    %Xticks{end} = Xlim_max;
    Xlabs = cellfun(@(x) sprintf('%+i',round(((x-Tzero)*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
    set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    set(gca,'XLim',[0 Xlim_max]);
    y_range = get(gca,'YLim');
    set(gca,'YLim',[0 max(y_range)]);
    set(gca,'Fontsize',fontsz)
    ylabel(Yax_labs{figidx})
    legend(simlabels)
    print(fullfile(fig_dir,[fig_fn '_onefig_' figIDs{figidx}]),'-djpeg')
end







function spiketrain = sim_spiketrain(sim_name,rtSpikes)

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
spiketrain = summed_timecourses(:,:,rtSpikes);


%stateswich_counts = cellfun(@(x) x{2}, sim_timecourse_output, 'UniformOutput', false);
%stateswich_counts = sum(cell2mat(stateswich_counts));

%spiketrain = summed_timecourses ./ stateswich_counts;

% %spikerate for 2ms bins for each cell
%num_binsamps = 2e-3/timestep; %num samples in 2ms
%cell_spkrates = spiketrain(:,:,rtSpikes);
%bin_magic = [numel(cell_spkrates(:,1)), num_binsamps, numel(cell_spkrates(1,:))/num_binsamps]; %set up for a magic trick
%cell_spkrates = reshape(cell_spkrates,bin_magic);
%cell_spkrates = squeeze(sum(cell_spkrates,2)) ./ (num_binsamps * timestep); %convert to Hz
%cell_spkrates = repmat(cell_spkrates,[num_binsamps 1 1]);
%spiketrain(:,:,rtSpikes) = reshape(cell_spkrates,size(spiketrain(:,:,rtSpikes))); %put the rabit back in the hat

end


