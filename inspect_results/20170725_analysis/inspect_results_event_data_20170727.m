clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.
%1) run combfig 20170725 prior to this if you're pulling from txtdump data.
%that'll make a datafile you can use here!
%2) this one seperates out the simulation data into different event windows
%to be analyzed by a seperate file. Also plots figures for the different
%events.
home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/';
addpath(fullfile(home_dir,'helper_functions'));
output_dir = fullfile(home_dir,'Results','change_switchspeed');
output_dir = fullfile(output_dir,'event_data');

if ~isdir(output_dir), mkdir(output_dir); end

how2load = 'txtdump'; %'txtdump' | 'regular' (will not covnert data from txtdump!!
sims2load = {'fastswitch_baseline','slowswitch_baseline',...
    'fastswitch_stimulus','slowswitch_stimulus'};

timestep = .25e-3;
recorded_switchtime =  250e-3; %actual switchtime in recorded switch
recorded_switchtime = recorded_switchtime/timestep; %actual switchtime in recorded switch

% %recorded vars (from noisycurrent_model()):
% %see vars2record code template comments at the very bottom of this file. it's basically legacy BS at this point
% num_vars2record = 1; %just have spikes this time
rtNoise = NaN;rtSpikes = 1; %just so I don't loose track of matrix inds
% var_inds = [rtSpikes]; %this is dumb just go with it, don't wana loose track
% fixvars = [4]; %only spikes  %fixvars = [1 4]; %spikes & noise

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);

num_sims = numel(sims2load);

spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    spike_data{loadidx} = sim_spikerate(sims2load{loadidx},timestep,rtSpikes,how2load);
end



num_events = 3; %prestim, stimulus preswitch, post-switch
event_labels = {'stable','pre-switch','post-switch'};
event_windows = cell(num_events,1); %event time bins
event_windows{1} = [-250e-3 , -150e-3]; % duration window, signage follows (T0 + X)
event_windows{2} = [-105e-3 , -5e-3];
event_windows{3} = [5e-3 , 105e-3];
event_windows = cellfun(@(x) x ./ timestep,event_windows,'UniformOutput',false); %convert to samples
event_windows = cellfun(@(x) (1+recorded_switchtime+x(1)):(recorded_switchtime+x(2)),...
    event_windows,'UniformOutput',false); %get window inds
adjusted_switch_onsets = cellfun(@(x) 1+recorded_switchtime-min(x),event_windows); %adjusted to the new plotting window

event_data = cell(size(event_windows));
for idx = 1:num_events
    curr_win = event_windows{idx};
    event_data{idx} = cellfun(@(x) x(:,curr_win),spike_data,'UniformOutput',false);
end


switch how2load
    case 'txtdump'
        %this indexing is from dump_data() function. kinda hardcoded but w/e
        %this also must match the indexing in make_datafile()
        cell_inds.Estay = 1;
        cell_inds.Istay = 2;
        cell_inds.Eswitch = 3;
        cell_inds.Iswitch = 4;
        %get rid of everything except spikerates
        for fix_idx = 1:num_events
            event_data{fix_idx} = cellfun(@(x) x(1:4,:), event_data{fix_idx},'UniformOutput',false);
        end
    case 'regular'
        %normal people indexing that makes sense
        cell_inds.Estay = celltype.excit & celltype.pool_stay;
        cell_inds.Eswitch = celltype.excit & celltype.pool_switch;
        cell_inds.Istay = celltype.inhib & celltype.pool_stay;
        cell_inds.Iswitch = celltype.inhib & celltype.pool_switch;
end


info.event_labels = event_labels;
info.sim_names = sims2load;
info.timestep = timestep;
info.recorded_switchtime = recorded_switchtime;
info.cell_inds = cell_inds;

save(fullfile(output_dir,'event_data'),'info','event_data','event_windows','adjusted_switch_onsets')


%plot for fun
fontsz = 12;
legend_labels = {'E-stay','E-switch','I-stay','I-switch'}; %ensure this lines up with how they're plotted
adj_legend = [-.25,-.055;-.25,0;0,.07];
simlabels = {'fastswitch baseline','slowswitch baseline',...
    'fastswitch stimulus','slowswitch stimulus   '};
Yax_labs = 'spike rate (Hz) per 2ms bin';

%event_labels = {'stable','pre-switch','post-switch'};
%.4 up a little more 

for figidx = 1:num_events
    
    data4plot = event_data{figidx};
    
    orient landscape
    for simidx = 1:num_sims
        
        subplot(num_sims/2,num_sims/2,simidx)
        
        %break it down by celltype & aggregate
        switch how2load
            case 'txtdump'
                %this indexing is from dump_data() function. kinda hardcoded but w/e
                %this also must match the indexing in make_datafile()
                Estay = 1;
                Istay = 2;
                Eswitch = 3;
                Iswitch = 4;
            case 'regular'
                %normal people indexing that makes sense
                Estay = celltype.excit & celltype.pool_stay;
                Eswitch = celltype.excit & celltype.pool_switch;
                Istay = celltype.inhib & celltype.pool_stay;
                Iswitch = celltype.inhib & celltype.pool_switch;
        end
        
        timecourse_data = data4plot{simidx};
        timecourse_data = timecourse_data(:,:,rtSpikes);
        
        Estay = mean(timecourse_data(Estay,:),1);
        Eswitch = mean(timecourse_data(Eswitch,:),1);
        Istay = mean(timecourse_data(Istay,:),1);
        Iswitch = mean(timecourse_data(Iswitch,:),1);
        
        
        %plot the aggregated timecourses
        hold on
        plot(Estay,'Linewidth',2)
        plot(Eswitch,'Linewidth',2)
        plot(Istay,'Linewidth',2)
        plot(Iswitch,'Linewidth',2)
        
        hold off
        Xlim_max = numel(timecourse_data(1,:));
        set(gca,'XLim',[0 Xlim_max]);
        Xticks = num2cell(get(gca,'Xtick'));
        if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
        adj_switch = adjusted_switch_onsets(figidx);
        Xlabs = cellfun(@(x) sprintf('%+i',((x-adj_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        if simidx == 3 || simidx == 4
            xlabel('state switch timecourse (ms)')
        end
        if simidx == 1 || simidx == 3
            ylabel(Yax_labs)
        end
        
        title(simlabels{simidx})
        
        %y_range = get(gca,'YLim');
        %x_range = get(gca,'XLim');
        %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
        set(gca,'Fontsize',fontsz)
        
        if simidx == 4
            leg = legend(legend_labels,'Fontsize',fontsz);
            legpos = get(leg,'Position'); %keep width & height
            legpos(1) = 1 - legpos(3);
            legpos(1) = legpos(1) + adj_legend(figidx,1);
            legpos(2) = legpos(2) + adj_legend(figidx,2);
            set(leg,'Position',legpos,'Units','normalized');
        end
        
    end
    fig_fn = strrep(event_labels{figidx},'-','');
    fig_fn = strrep(fig_fn,' ','_');
    print(fullfile(output_dir,[fig_fn ' _timecourse']),'-djpeg')
    close all
end




function avg_timecourse = sim_spikerate(sim_name,timestep,rtSpikes,how2load)

config_options.modeltype = 'JK';
config_options.sim_name = sim_name;
options = set_options(config_options);

%get results
switch how2load
    case 'regular'
        output_fn = [options.modeltype '_' options.sim_name '.mat'];
    case 'txtdump'
        output_fn = [options.modeltype '_' options.sim_name '_txtdump.mat'];
end
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


%vars2record code template!!
%
% num_vars2record = 2;
% rtNoise = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
% var_inds = [rtNoise,rtSpikes]; %this is dumb just go with it, don't wana loose track
%lets record, noise, Sg, D, spikes, (I think that's it?)
%num_vars2record = 4;
%rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
%var_inds = [rtNoise,rtSg,rtD,rtSpikes]; %this is dumb just go with it, don't wana loose track
%get all combinations of cross-correlations, 4 per simulation?