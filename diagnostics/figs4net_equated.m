clear
clc
format compact
hold off;close all

%show spikerate plots when slow & fast networks have equated  sampling time
%this requires saved datafiles from a "diagnostic" model run


addpath('../')

%specify what results
Sname = 'example_behavior_equated';
jobs = sort([113:10:193,114:10:194]); %do a few runs for net #2

%figure options
lnsz = 2;
fontsz = 20;
shift_start = 12; %plot starting x seconds into simulation data
Tmax_plot = 60; %0, or specify total plot duration (plot ends at Tmax_plot + shift_start)
window_sz = 150e-3; %for averaging spikerates

for j = 1:numel(jobs)
    
    close all
    
    options = set_options('modeltype','diagnostics','comp_location','woodstock',...
        'sim_name',Sname,'jobID',jobs(j));
        
    
    %cell_order = {'E-stay','E-switch','I-stay','I-switch'}; %plot these in order
    cell_order = {'E-stay','E-switch'}; %only plot excitatory
    Npools = numel(cell_order);
    cell_fns = strrep(cell_order,'-',''); %for structure fieldnames
    cell_cols = get(gca,'colororder');close all
    cell_cols = cell_cols(1:Npools,:);
    cell_cols = num2cell(cell_cols,2);
    
    %make a figure directory
    fig_dir = fullfile(options.save_dir,options.sim_name);
    if ~isdir(fig_dir),mkdir(fig_dir);end
    
    
    %create this stuff, not saved with data
    pool_options.num_cells = 250;
    pool_options.sz_pools = [.5 .5]; %proportion stay & switch
    pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
    pool_options.p_conn = .5; %connection probability 50%
    celltype = celltype_logicals(pool_options);
    
    %load the datafiles
    savename = fullfile(options.save_dir,options.sim_name);
    load(savename)
    load(sprintf('%s_D',savename))
    load(sprintf('%s_V',savename))
    load(sprintf('%s_spikes',savename))
    load(sprintf('%s_S',savename))
    timestep = options.timestep; %easier var name
    
    %get average state duration.. this won't usually capture last state but it's close enough
    durations = sim_results{1};
    if ~isempty(durations)
        %just get all of them, baseline test. Everything that's not undecided
        valid_states = startsWith(durations(:,end),'stim');
        if sum(valid_states) > 0
            durations = cell2mat(durations(valid_states,1:2));
            durations = durations * timestep;
            timecourse = array2table(durations,'VariableNames',{'event_time','duration'});
            timecourse.start = timecourse.event_time - timecourse.duration;
            timecourse = timecourse(timecourse.event_time > shift_start,:);
            timecourse = timecourse(timecourse.start > shift_start,:);
            if shift_start > 0
                trunc = timecourse.start < shift_start;
                timecourse.duration(trunc) = timecourse.duration(trunc) - ...
                    (shift_start - timecourse.start(trunc));
            end
            
            if Tmax_plot > 0
                timecourse = timecourse(timecourse.start <= Tmax_plot,:);
                trunc = timecourse.event_time > Tmax_plot;
                timecourse.duration(trunc) = timecourse.duration(trunc) - ...
                    (timecourse.event_time(trunc) - Tmax_plot);
            end
            fprintf('------------------------\n')
            disp(timecourse)
            fprintf('\nmean state duration = %.2fs\n\n',mean(timecourse.duration))
        end
    end
    
    if shift_start > 0
        start_ind = round(shift_start / timestep);
        %now truncate
        spikes = spikes(:,start_ind:end);
        Drec_fast = Drec_fast(:,start_ind:end);
        Drec_slow = Drec_slow(:,start_ind:end);
        Vrec = Vrec(:,start_ind:end);
        Srec = Srec(:,start_ind:end);
    end
    if Tmax_plot > 0
        Tmax_ind = 1:round(Tmax_plot / timestep);
        %now truncate
        spikes = spikes(:,Tmax_ind);
        Drec_fast = Drec_fast(:,Tmax_ind);
        Drec_slow = Drec_slow(:,Tmax_ind);
        Vrec = Vrec(:,Tmax_ind);
        Srec = Srec(:,Tmax_ind);
    end
    
    S = sim_windowrate(spikes,timestep,celltype,window_sz);
    
    %only need spikerate plots
    %spikerates
    hold on
    for idx = 1:Npools
        plot(S.(cell_fns{idx}),'Linewidth',lnsz,'Color',cell_cols{idx})
    end
    Nobs = unique(structfun(@numel,S));
    xlim([1,Nobs])
    Xticks = num2cell(get(gca,'Xtick'));
    Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
    Xlabs = strrep(Xlabs,'.0','');
    set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
    ylabel({'pool spiking (Hz)'},'Fontweight','b')
    xlabel('seconds','Fontweight','b');
    set(gca,'box','off');
    set(gca,'FontSize',fontsz);
    % hold on %turn off box to remove upper/right ticks & add border back in
    % border = linspace(min(get(gca,'xlim')),max(get(gca,'xlim')));
    % plot(border,zeros(size(border))+max(get(gca,'ylim')),'color','k')
    % border = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')));
    % plot(zeros(size(border))+max(get(gca,'xlim')),border,'color','k')
    axP = get(gca,'Position');
    [lp, ~] = legend(cell_order,'FontWeight','b',...
        'Location','northoutside','Box','off','Orientation','horizontal','FontSize',fontsz);
    set(gca, 'Position', axP)
    lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];
    
    fprintf('printing figure for %s\n',options.sim_name)
    FN = fullfile(fig_dir,'spikerates_Eonly.png');
    print(FN,'-dpng','-r600');
    savefig(gcf,strrep(FN,'.png',''),'compact')
end

function cell_data = sim_spikerate(cell_raster,timestep,celltype)

binsz = 75e-3;
num_binsamps = binsz./timestep; %num samples in Xms (the bin size)
if round(num_binsamps) - num_binsamps < 1e-12 %roundoff errors...
    num_binsamps = round(num_binsamps);
end
if mod(num_binsamps,1) ~= 0 %not evenly divisible... find next best thing
    alt_binsz = binsz-(5e-3):1e-3:binsz+(5e-3);
    even_div = rem(alt_binsz, timestep) == 0;
    alt_binsz = alt_binsz(even_div);
    [~,binsz] = min(abs(binsz - alt_binsz)); %find closest evenly divisible binsize
    binsz = alt_binsz(binsz);
    num_binsamps = binsz/timestep;
end

raster_sz = size(cell_raster);
if mod(raster_sz(2),num_binsamps) ~= 0 %you have to trim it down, equally divisible by bin size
    cell_raster = cell_raster(:,1:end - mod(raster_sz(2),num_binsamps));
    raster_sz = size(cell_raster); %should be good now
end
bin_magic = [raster_sz(1), num_binsamps,raster_sz(2)/num_binsamps]; %set up for a magic trick
cell_raster = reshape(cell_raster,bin_magic);
cell_raster = squeeze(sum(cell_raster,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_raster = repmat(cell_raster,[num_binsamps 1 1]);
cell_raster = reshape(cell_raster,raster_sz); %put the rabit back in the hat

%normal people indexing that makes sense, then take the mean
cell_data.Estay = mean(cell_raster(celltype.excit & celltype.pool_stay,:),1);
cell_data.Eswitch = mean(cell_raster(celltype.excit & celltype.pool_switch,:),1);
cell_data.Istay = mean(cell_raster(celltype.inhib & celltype.pool_stay,:),1);
cell_data.Iswitch = mean(cell_raster(celltype.inhib & celltype.pool_switch,:),1);
end

function cell_data = sim_windowrate(cell_raster,timestep,celltype,window_sz)

num_binsamps = window_sz/timestep; %num samples in Xms
if round(num_binsamps) - num_binsamps < 1e-12 %roundoff errors...
    num_binsamps = round(num_binsamps);
end
if mod(num_binsamps,1) ~= 0 %not evenly divisible... find next best thing
    alt_binsz = window_sz-(5e-3):1e-3:window_sz+(5e-3);
    even_div = rem(alt_binsz, timestep) == 0;
    alt_binsz = alt_binsz(even_div);
    [~,binsz] = min(abs(binsz - alt_binsz)); %find closest evenly divisible binsize
    binsz = alt_binsz(binsz);
    num_binsamps = binsz/timestep;
end

cell_raster = num2cell(cell_raster,2);
k = ones(1, num_binsamps);
k = k ./ (num_binsamps * timestep); %convert to Hz
cell_raster = cellfun(@(x) conv(x, k, 'same'),cell_raster,'UniformOutput',false);
cell_raster = cat(1,cell_raster{:});

%normal people indexing that makes sense, then take the mean
cell_data.Estay = mean(cell_raster(celltype.excit & celltype.pool_stay,:),1);
cell_data.Eswitch = mean(cell_raster(celltype.excit & celltype.pool_switch,:),1);
cell_data.Istay = mean(cell_raster(celltype.inhib & celltype.pool_stay,:),1);
cell_data.Iswitch = mean(cell_raster(celltype.inhib & celltype.pool_switch,:),1);
end

function cell_data = simple_pool_avg(x,celltype)
%takes data matrix and returns simple pool averages
%use for naturally smoothed variables, depression, syn gating, etc
cell_data.Estay = mean(x(celltype.excit & celltype.pool_stay,:),1);
cell_data.Eswitch = mean(x(celltype.excit & celltype.pool_switch,:),1);
cell_data.Istay = mean(x(celltype.inhib & celltype.pool_stay,:),1);
cell_data.Iswitch = mean(x(celltype.inhib & celltype.pool_switch,:),1);
end


function [on_bounds,off_bounds] = get_state_limits(S)
%takes the spikerate structre S
%returns the Y limits for plotting the active & inactive pool activity
I = [S.Istay;S.Iswitch]; E = [S.Estay;S.Eswitch];

[Emin,Emax] = bounds(E); [Imin,Imax] = bounds(I);

maxlim = @(x,y) max(mean(x),mean(y)) + 3*max(std(x),std(y));
minlim = @(x,y) min(mean(x),mean(y)) - 3*min(std(x),std(y));

off_bounds(1) = 0; %lets default to this...
off_bounds(2) = 9; %off_bounds(2) = maxlim(Emin,Imin);

on_bounds(1) = minlim(Emax,Imax);
on_bounds(2) = maxlim(Emax,Imax);
end

function r = get_raster_img(S,dt,Sdur)
%takes spike matrix, timestep, and
%how wide you want spikes shown (duration, in seconds). Outputs raster

kernel = (Sdur / dt);
if mod(floor(kernel),2) == 1
    kernel = floor(kernel);
elseif mod(ceil(kernel),2) == 1
    kernel = ceil(kernel);
else
    kernel = kernel + 1; %must be an even integer
end
kernel = ones(1,kernel);


r = NaN(size(S));
for idx = 1:numel(S(:,1)) %for each cell
    r(idx,:) = conv(S(idx,:),kernel,'same');
end
r(r > 1) = 1; %just in case spikes were right next to each other or something


end

