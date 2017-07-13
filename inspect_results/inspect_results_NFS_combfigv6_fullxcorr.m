clear
clc
close all
format compact

%NOTE: this script requires 2017a. Local function below.

fig_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/Results/NFSv5_combfigs2';
fig_dir = fullfile(fig_dir,'fullXcorr_long_coeff');
if ~isdir(fig_dir)
    mkdir(fig_dir)
end

how2load = 'txtdump'; %'txtdump' | 'regular'
min_time = 1.25; %min switch time in seconds (for txtdump function)
max_time = 9999; %max switch time in seconds (for txtdump function)
name4plots = 'late_timecourse'; %added to basic fig directory, distingish different timecourse plots etc
recorded_switchtime =  250e-3; %actual switchtime in recorded switch
preswitch_plottime = 105e-3; %preswitch duration to plot (T0-X)
postswitch_plottime = -5e-3; %postwitch duration to plot (T+X)

scaled_plot = 'off'; %shouldn't need this anymore
scalefac = 2;
%result summaries
fontsz = 12;
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..
%recorded vars (from noisycurrent_model()):
num_vars2record = 1; %just have spikes this time
rtNoise = NaN;rtSpikes = 1; %just so I don't loose track of matrix inds
var_inds = [rtSpikes]; %this is dumb just go with it, don't wana loose track
fixvars = [4]; %only spikes  %fixvars = [1 4]; %spikes & noise
%these should really all line up... pay attn to that
sims2load = {'NFSv5_Estay','NFSv5_Eswitch','NFSv5_Istay','NFSv5_Iswitch'};
simlabels = {'Exit. stay (-)','Exit. switch (+)','Inhib. stay (+)','Inhib. switch (-)'};
legend_labels = {'E-stay','E-switch','I-stay','I-switch'};
% num_vars2record = 2;
% rtNoise = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
% var_inds = [rtNoise,rtSpikes]; %this is dumb just go with it, don't wana loose track
%lets record, noise, Sg, D, spikes, (I think that's it?)
%num_vars2record = 4;
%rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
%var_inds = [rtNoise,rtSg,rtD,rtSpikes]; %this is dumb just go with it, don't wana loose track
%get all combinations of cross-correlations, 4 per simulation?
num_sims = numel(sims2load);

% switch how2load
%     case 'txtdump'
%         for loadidx = 1:num_sims
%             make_datafile(sims2load{loadidx},timestep,min_time,max_time);
%         end
% end

spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    spike_data{loadidx} = sim_spikerate(sims2load{loadidx},timestep,rtSpikes,how2load);
end

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%set up for muiltiple figs
fig_fn = 'switching_Xcorr_full';
figIDs = {'noise','Sg','D','spikes'};
%base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','pool-1 * pool-2'};
%I only did noise & spikes, just index these so I dont have to copy & paste etc..
figIDs = figIDs(fixvars);
legloc = legloc(fixvars);
Yax_labs = Yax_labs(fixvars);


%only plot -Xms to +Xms
recorded_switchtime = recorded_switchtime/timestep; %actual switchtime in recorded switch
postswitch_plottime = postswitch_plottime/timestep;
preswitch_plottime = preswitch_plottime/timestep;
record_duration = size(cat(1,spike_data{:}),2); %get the duration of recorded timecourses
plotting_window = 1 + recorded_switchtime - preswitch_plottime:recorded_switchtime + postswitch_plottime;
onset_switch = 1 + recorded_switchtime - min(plotting_window); %adjusted to the new plotting window
%cut down the data matrix to this window
spike_data = cellfun(@(x) x(:,plotting_window,:),spike_data,'UniformOutput',false);

%get all combinations of cross-correlations, 4 per simulation?
switch how2load
    case 'txtdump'
        %this indexing is from dump_data() function. kinda hardcoded but w/e
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

celltype_inds = {Estay,Eswitch,Istay,Iswitch};
celltype_labels = {'E-stay','E-switch','I-stay','I-switch'};
Xcor_pairings = nchoosek(1:4,2);
Xcor_pairings = num2cell(Xcor_pairings,2);
num_Xcors = numel(Xcor_pairings);
Xcor_inds = cellfun(@(x) celltype_inds(x),Xcor_pairings,'UniformOutput',false);
Xcor_labels = cellfun(@(x) celltype_labels(x),Xcor_pairings,'UniformOutput',false);

Rcells = cell(num_Xcors,num_sims);
for simidx = 1:num_sims
    sim_data = spike_data{simidx};
    for XCidx = 1:num_Xcors
        curr_pair = Xcor_inds{XCidx};
        curr_data = cellfun(@(x) sim_data(x,:),curr_pair,'UniformOutput',false);
        curr_data = cellfun(@(x) mean(x,1),curr_data,'UniformOutput',false); %mean spikerate across cells
        curr_data = cellfun(@(x) x - mean(x),curr_data,'UniformOutput',false); %zero mean data   
        %[Rcells{XCidx,simidx},laginds] = xcorr(curr_data{1},curr_data{2},'unbiased'); %do Xcorr
        [Rcells{XCidx,simidx},laginds] = xcorr(curr_data{1},curr_data{2},'coeff'); %do Xcorr
    end
end


for figidx = 1:num_vars2record
    
    %subplot(num_sims/2,num_sims/2,num_sims)
    figure(1)
    orient landscape
    for simidx = 1:num_sims
        
        subplot(num_sims/2,num_sims/2,simidx)
        
        for XCidx = 1:num_Xcors
            hold on
            plot(laginds,Rcells{XCidx,simidx},'Linewidth',2)
        end
        hold off
        set(gca,'XLim',[min(laginds)-1 max(laginds)+1]);
        Xticks = num2cell(get(gca,'Xtick'));
        Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        
        if simidx == 3 || simidx == 4
            xlabel('pool-2 lag offset (ms)')
        end
        if simidx == 1 || simidx == 3
            ylabel(Yax_labs{figidx})
        end
        %         if simidx == 4
        %             leg_labels = cellfun(@(x) [x{1} ' & ' x{2}],Xcor_labels,'UniformOutput', false);
        %             leg = legend(leg_labels,'Fontsize',fontsz);
        %             legpos = get(leg,'Position'); %keep width & height
        %             legpos(1) = 1 - legpos(3);
        %             legpos(2) = legpos(2) + .07;
        %             set(leg,'Position',legpos,'Units','normalized');
        %         end
        title(simlabels{simidx})
        set(gca,'Fontsize',fontsz)
    end
    
    %align all y axes
    allax = gca;
    allax = allax.Parent.Children;
    linkaxes(allax,'y');
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx}]),'-djpeg')
    
    
    %print a stupid external legend
    close all;hold off
    figure(1)
    plot(repmat(1:num_Xcors,num_Xcors,1),'linewidth',2);
    leg_labels = cellfun(@(x) [x{1} ' & ' x{2}],Xcor_labels,'UniformOutput', false);
    leg = legend(leg_labels,'Fontsize',fontsz+10);
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx} '_legend']),'-djpeg') %stupid legend
    
    
    
    close all;hold off
    figure(2)
    orient tall
    fig2_titles = cellfun(@(x) [x{1} ' & ' x{2}],Xcor_labels,'UniformOutput', false);
    for XCidx = 1:num_Xcors
        subplot(num_Xcors/2,2,XCidx)
        for simidx = 1:num_sims
            
            hold on
            plot(laginds,Rcells{XCidx,simidx},'Linewidth',2)
        end
        hold off
        set(gca,'XLim',[min(laginds)-1 max(laginds)+1]);
        Xticks = num2cell(get(gca,'Xtick'));
        Xlabs = cellfun(@(x) sprintf('%+i',round((x*timestep)/1e-3)),Xticks,'UniformOutput', false); %this is for normal stuff
        set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
        if XCidx == 5 || XCidx == 6
            xlabel('pool-2 lag offset (ms)')
        end
        if mod(XCidx,2) == 1
            ylabel(Yax_labs{figidx})
        end
        title(fig2_titles{XCidx})
        set(gca,'Fontsize',fontsz)
    end
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx} '2']),'-djpeg')

    %print a stupid external legend
    close all;hold off
    figure(1)
    plot(repmat(1:num_sims,num_sims,1),'linewidth',2);
    leg = legend(simlabels,'Fontsize',fontsz+10);
    print(fullfile(fig_dir,[fig_fn '_' figIDs{figidx} '2_legend']),'-djpeg') %stupid legend
    close all;hold off
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


function make_datafile(sim_name,timestep,Tmin,Tmax)

config_options.modeltype = 'JK';
config_options.sim_name = sim_name;
options = set_options(config_options);
dump_dir = fullfile(options.save_dir,'timecourse_dump');
fns = dir(fullfile(dump_dir,'*txt'));
fns = {fns.name}';
%get switch durations from filenames
switch_durs = cellfun(@(x) strsplit(x,'_'),fns,'UniformOutput',false);
switch_durs = cellfun(@(x) x{end},switch_durs,'UniformOutput',false);
switch_durs = cellfun(@(x) strrep(x,'.txt',''),switch_durs,'UniformOutput',false);
switch_durs = cellfun(@(x) str2double(x),switch_durs,'UniformOutput',false);
switch_durs = cell2mat(switch_durs);
%find which files match input criteria
Tmax = Tmax / timestep; %convert to samples
Tmin = Tmin / timestep;

valid_data = switch_durs <= Tmax & switch_durs >= Tmin;
num_switches = sum(valid_data);
valid_files = fns(valid_data);

data = load(fullfile(dump_dir,valid_files{1})); %get dims
data = zeros(size(data)); %sum each file, too much for memory...
for fnidx = 1:num_switches
    if mod(fnidx,1000) == 0,fprintf('loading file #%i/%i...\n',fnidx,num_switches);end
    curr_file = load(fullfile(dump_dir,valid_files{fnidx}));
    data = data + curr_file;
end

sim_switch_timecourses = {[{data},{num_switches}]}; %mimic simulation output file
output_fn = [options.modeltype '_' options.sim_name '_txtdump.mat'];
save(fullfile(options.save_dir,output_fn),'sim_switch_timecourses')
end


function muC = zeromean(X)
muC = X - mean(X);
end

