clear
clc
close all
format compact

%MORE IMPORTANT NOTE: this code looks like it was messed with... 
%NOTE: this script requires 2017a. Local function below. 

scaled_plot = 'off';
scalefac = 2;
fig_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/Results/NFSv4_combfigs';
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

sims2load = {'NFSv4_Estay','NFSv4_Eswitch','NFSv4_Iswitch','NFSv4_Istay'};
simlabels = {'Exit. stay (-)','Exit. switch (+)','Inhib. switch (-)','Inhib. stay (+)'};
num_sims = numel(sims2load);

spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    curr_data = sim_spikerate(sims2load{loadidx},timestep,rtSpikes);
    spike_data{loadidx} = curr_data(:,:,rtSpikes);
end

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 150;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [2/3 1/3]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%set up for muiltiple figs
fig_fn = 'switching_timecourse';
figIDs = {'noise','Sg','D','spikes'};
%base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','spike rate (Hz) per 2ms bin'};
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
    
    
    %subplot(num_sims/2,num_sims/2,num_sims)
    figure(1)
    orient landscape
    for simidx = 1:num_sims
        
        subplot(num_sims/2,num_sims/2,simidx)
        
        %break it down by celltype & aggregate
        Estay = celltype.excit & celltype.pool_stay;
        Eswitch = celltype.excit & celltype.pool_switch;
        Istay = celltype.inhib & celltype.pool_stay;
        Iswitch = celltype.inhib & celltype.pool_switch;
        
        timecourse_data = spike_data{simidx};
        
        Estay = mean(timecourse_data(Estay,:));
        Eswitch = mean(timecourse_data(Eswitch,:));
        Istay = mean(timecourse_data(Istay,:));
        Iswitch = mean(timecourse_data(Iswitch,:));
        
        switch scaled_plot
            case 'off'
                %plot the aggregated timecourses
                hold on
                plot(Estay,'Linewidth',2)
                plot(Eswitch,'Linewidth',2)
                plot(Istay,'Linewidth',2)
                plot(Iswitch,'Linewidth',2)
                hold off
                
                Xlim_max = numel(timecourse_data(1,:)); %this is really annoying, probably remove this if timecourse changes and it's unneeded. 
                onset_switch = 250e-3/timestep; %add the switching time
                Xticks = num2cell(get(gca,'Xtick'));
                Xticks{end} = Xlim_max;
                Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
                set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
               
                if figidx == rtNoise %just do for noise..
                    Yticks = num2cell(get(gca,'Ytick'));
                    Ylabs = cellfun(@(x) sprintf('%.0f',x/1e-12),Yticks,'UniformOutput', false);
                    set(gca, 'YTickLabel', Ylabs,'Ytick',cell2mat(Yticks));
                end
                if simidx == 3 || simidx == 4
                    xlabel('state switch timecourse (ms)')
                end
                ylabel(Yax_labs{figidx})
                title(simlabels{simidx})
                
                %y_range = get(gca,'YLim');
                %x_range = get(gca,'XLim');
                set(gca,'XLim',[0 Xlim_max]);
                %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
                set(gca,'Fontsize',fontsz)
                
                
                
                %(only for spikerates)
                %plot the spikerate timecourses scaled down by some factor
            case 'on'
                hold on
                plot(Estay ./ scalefac,'Linewidth',2)
                plot(Eswitch ./ scalefac,'Linewidth',2)
                plot(Istay ./ scalefac,'Linewidth',2)
                plot(Iswitch ./ scalefac,'Linewidth',2)
                hold off
                
                Xlim_max = numel(timecourse_data(1,:)); %this is really annoying, probably remove this if timecourse changes and it's unneeded.
                onset_switch = 250e-3/timestep; %add the switching time
                Xticks = num2cell(get(gca,'Xtick'));
                Xticks{end} = Xlim_max;
                Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
                
                set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
                if simidx == 3 || simidx == 4
                    xlabel('state switch timecourse (ms)')
                end
                ylabel(Yax_labs{figidx})
                title(simlabels{simidx})
                %y_range = get(gca,'YLim');
                %x_range = get(gca,'XLim');
                set(gca,'XLim',[0 Xlim_max]);
                %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
                set(gca,'Fontsize',fontsz)
                
        end
        
        if simidx == 3
            legend({'E-stay','E-switch','I-stay','I-switch'},'location','northwest',...
                'Fontsize',fontsz)
        end
    end
    switch scaled_plot
        case 'off'
            print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx}]),'-djpeg')
        case 'on'
            print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx} '_scaled']),'-djpeg')
    end
end

function avg_timecourse = sim_spikerate(sim_name,timestep,rtSpikes)
% MYMEDIAN Another example of a local function.

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


