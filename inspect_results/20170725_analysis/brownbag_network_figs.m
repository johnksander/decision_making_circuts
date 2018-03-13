clear
clc
close all
format compact

%this comes from writeup_inspect_results_combfig_20170725.m
%changing the figures a bit for brownbag
%NOTE: this script requires 2017a. Local function below.
home_dir = '/Users/ksander/Desktop/work/ACClab/rotation/project/';
addpath(home_dir)
addpath(fullfile(home_dir,'helper_functions'))
fig_dir = fullfile(home_dir,'Results/brownbag');
if ~isdir(fig_dir),mkdir(fig_dir);end
fig_name = 'combfig';

conv_txtdump = 'off'; %off/on
how2load = 'txtdump'; %'txtdump' | 'regular'
min_time = .6; %min switch time in seconds (for txtdump function)
max_time = 9999; %max switch time in seconds (for txtdump function)
name4plots = 'combfig'; %added to basic fig directory, distingish different timecourse plots etc
recorded_switchtime =  250e-3; %actual switchtime in recorded switch
%preswitch_plottime = 250e-3; %preswitch duration to plot (T0-X)
%postswitch_plottime = 150e-3; %postwitch duration to plot (T+X)
preswitch_plottime = 105e-3; %preswitch duration to plot (T0-X)
postswitch_plottime = -5e-3; %postwitch duration to plot (T+X)


%scaled_plot = 'off'; %shouldn't need this anymore
scalefac = 2;
%result summaries
fontsz = 26;
stim_labels = {'stim A','stim B'};
timestep = .25e-3; %this should really make it's way into set_options(), used for conv2secs here..
%recorded vars (from noisycurrent_model()):
num_vars2record = 1; %just have spikes this time
rtNoise = NaN;rtSpikes = 1; %just so I don't loose track of matrix inds
var_inds = [rtSpikes]; %this is dumb just go with it, don't wana loose track
fixvars = [4]; %only spikes  %fixvars = [1 4]; %spikes & noise
%these should really all line up... pay attn to that
%sims2load = {'fastswitch_baseline','slowswitch_baseline',...
%    'fastswitch_stimulus','slowswitch_stimulus'};
sims2load = {'fastswitch_stimulus','slowswitch_stimulus'};
simlabels = {'Entice to stay','Repel to leave'};
legend_labels = {'E-stay','E-switch','I-stay','I-switch'};

% num_vars2record = 2;
% rtNoise = 1;rtSpikes = 2; %just so I don't loose track of matrix inds
% var_inds = [rtNoise,rtSpikes]; %this is dumb just go with it, don't wana loose track
%lets record, noise, Sg, D, spikes, (I think that's it?)
%num_vars2record = 4;
%rtNoise = 1;rtSg = 2;rtD = 3;rtSpikes = 4; %just so I don't loose track of matrix inds
%var_inds = [rtNoise,rtSg,rtD,rtSpikes]; %this is dumb just go with it, don't wana loose track

orange = [250 70 22]./255;
matblue = [0,0.4470,0.7410];
fig_dir = fullfile(fig_dir,name4plots);
if ~isdir(fig_dir),mkdir(fig_dir),end

%remake celltype logicals.. (if you use this code again, check this over!!!!!)
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


num_sims = numel(sims2load);

switch conv_txtdump
    case 'on'
        for loadidx = 1:num_sims
            make_datafile(sims2load{loadidx},min_time,max_time,timestep,celltype);
        end
end


spike_data = cell(num_sims,1);
for loadidx = 1:num_sims
    spike_data{loadidx} = sim_spikerate(sims2load{loadidx},timestep,rtSpikes,how2load);
end


%set up for muiltiple figs
%fig_fn = 'switching_timecourse_%i';
fig_fn = 'preswitch_timecourse_%i';
figIDs = {'noise','Sg','D','spikes'};
%base_title = {'cell-type mean timecourse';'stimuli A (1.5 x current) to switch-state (1 x current)'};
%legloc = {'southwest','northwest','West','southwest'};
legloc = {'west','northwest','West','west'};
%Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','spike rate (Hz) per 2ms bin'};
Yax_labs = {'current noise (picoamps)','synaptic gating','synaptic depression','spike rate (Hz)'};
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


figidx = 1;%:num_vars2record
close all
hold off

orient landscape
for simidx = 1:num_sims %this is the actual figidx now.. not doing subplots
    figure
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
    
    timecourse_data = spike_data{simidx};
    timecourse_data = timecourse_data(:,:,var_inds(figidx));
    
    Estay = mean(timecourse_data(Estay,:),1);
    Eswitch = mean(timecourse_data(Eswitch,:),1);
    Istay = mean(timecourse_data(Istay,:),1);
    Iswitch = mean(timecourse_data(Iswitch,:),1);
    
    %plot the aggregated timecourses
    hold on
    lnsz = 5;
    plot(Estay,'Linewidth',lnsz)
    plot(Eswitch,'Linewidth',lnsz)
    plot(Istay,'Linewidth',lnsz)
    plot(Iswitch,'Linewidth',lnsz)
    
    hold off
    Xlim_max = numel(timecourse_data(1,:));
    set(gca,'XLim',[0 Xlim_max]);
    Xticks = num2cell(get(gca,'Xtick'));
    
    
    %Xticks = Xticks(2:2:end-1);
    if numel(Xticks) > 5,Xticks = Xticks(1:2:numel(Xticks)); end %for crowded axes
    Xlabs = cellfun(@(x) sprintf('%+i',((x-onset_switch)*timestep)/1e-3),Xticks,'UniformOutput', false); %this is for normal stuff
    set(gca, 'XTickLabel', Xlabs,'Xtick',cell2mat(Xticks));
    
    xlabel('Leave decision (ms)','FontWeight','b')
    ylabel(Yax_labs{figidx},'FontWeight','b')
    title(simlabels{simidx})
    
    %y_range = get(gca,'YLim');
    %x_range = get(gca,'XLim');
    %text(.025*max(x_range),.9*max(y_range),sprintf('n switches = %i',stateswich_counts))
    set(gca,'Fontsize',fontsz)
    
    
    %(only for spikerates)
    %plot the spikerate timecourses scaled down by some factor
    
    
    %     if simidx == 4
    %         leg = legend(legend_labels,'Fontsize',fontsz);
    %         legpos = get(leg,'Position'); %keep width & height
    %         legpos(1) = 1 - legpos(3);
    %         legpos(2) = legpos(2) + .07;
    %         set(leg,'Position',legpos,'Units','normalized');
    %     end
    
    
    print(fullfile(fig_dir,sprintf(fig_fn,simidx)),'-djpeg')
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

function make_datafile(sim_name,Tmin,Tmax,timestep,celltype)

fprintf('----loading textdump data for simulation: %s\n',sim_name)

config_options.modeltype = 'JK';
config_options.sim_name = sim_name;
options = set_options(config_options);
dump_dir = fullfile(options.save_dir,'timecourse_dump');
%just do both stims for now!!
fns = dir(fullfile(dump_dir,'*','*txt'));
%fns = {fns.name}';
fns = cellfun(@(x,y) fullfile(x,y),{fns.folder}',{fns.name}','UniformOutput',false);
%get switch durations from filenames
switch_durs = cellfun(@(x) strsplit(x,'.txt'),fns,'UniformOutput',false);
switch_durs = cellfun(@(x) x{1},switch_durs,'UniformOutput',false);
switch_durs = cellfun(@(x) strsplit(x,'_'),switch_durs,'UniformOutput',false);
switch_durs = cellfun(@(x) str2double(x{end}),switch_durs,'UniformOutput',false);
switch_durs = cat(1,switch_durs{:});
%find which files match input criteria
Tmax = Tmax / timestep; %convert to samples
Tmin = Tmin / timestep;

valid_data = switch_durs <= Tmax & switch_durs >= Tmin;
num_switches = sum(valid_data);
valid_files = fns(valid_data);

fprintf('----valid switches found = %i\n',num_switches)

data = load(valid_files{1}); %get dims
data = zeros(size(data)); %sum each file, too much for memory...

for fnidx = 1:num_switches
    if mod(fnidx,1000) == 0,fprintf('loading file #%i/%i...\n',fnidx,num_switches);end
    curr_file = load(valid_files{fnidx});
    data = data + curr_file;
end

%need to convert to from aggragate spikes to average cell spikes or
%everything gets messed up downstream...

%this indexing is from dump_data() function. kinda hardcoded but w/e
%must match indexing outside this function!!
Estay = 1;
Istay = 2;
Eswitch = 3;
Iswitch = 4;

%this is also super ugly code right here... again, w/e just roll with it
num_Estay = sum(celltype.excit & celltype.pool_stay);
num_Eswitch = sum(celltype.excit & celltype.pool_switch);
num_Istay = sum(celltype.inhib & celltype.pool_stay);
num_Iswitch = sum(celltype.inhib & celltype.pool_switch);

Ctype_inds = {Estay,Eswitch,Istay,Iswitch};
num_cells = {num_Estay,num_Eswitch,num_Istay,num_Iswitch};
for fix_idx = 1:numel(Ctype_inds)
    data(Ctype_inds{fix_idx},:) = data(Ctype_inds{fix_idx},:) ./ num_cells{fix_idx};
end

%finally, save and try not to think about the messy code above

sim_switch_timecourses = {[{data},{num_switches}]}; %mimic simulation output file
output_fn = [options.modeltype '_' options.sim_name '_txtdump.mat'];
save(fullfile(options.save_dir,output_fn),'sim_switch_timecourses')

end






