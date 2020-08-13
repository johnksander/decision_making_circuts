clear
clc
format compact
hold off;close all

%investigating model behavior

addpath('../')
jobs = 12:14;
sname = 'test_nodata';
do_mode = 'load'; %'run' | 'load'

Njobs = numel(jobs);
data = cell(Njobs,1);

for idx = 1:numel(jobs)
    %my model
    %---setup---------------------
    options = set_options('modeltype','diagnostics','comp_location','bender',...
        'sim_name','test_nodata','jobID',jobs(idx));
    
    switch do_mode
        case 'run'
            fprintf('loading job JID = %i...',jobs(idx))
            data{idx} = get_data(options);
            fprintf('finished\n')
    end
    
end


fig_dir = fullfile(options.save_dir,'special_inspect');
if ~isdir(fig_dir),mkdir(fig_dir);end
FN = fullfile(fig_dir,'data.mat');

switch do_mode
    case 'run'
        save(FN,'data')
    case 'load'
        load(FN)
end

lnsz = 2; %spikerate plots
fontsz = 12;

% mean synaptic output variable, mean firing rate, mean depression variable, 
%the product of rate x depression, etc., all as a function of time since stimulus 
%onset for the Stay pool in the case of the standard (purely hedonic) input then in 
%the case of increased input.

do_cells = {'E-stay'};
Npools = numel(do_cells);
opts = cellfun(@(x) x.options,data,'UniformOutput',false);
timestep = opts{1}.timestep;
stims = cellfun(@(x) x.trial_stimuli{1} ,opts,'UniformOutput',false);
%this is a little hardcoded!
stims = cellfun(@unique,stims);

%reorder the data based on increasing stim intensity 
[stims,new_order] = sort(stims);
opts = opts(new_order);
data = data(new_order);

leg_labs = cellfun(@(x) sprintf('%.0f Hz Estay',x),num2cell(stims),'UniformOutput',false);

%1) synaptic output
Ylab = {'E-stay','synaptic gating'};
figure
for idx = 1:Njobs
    curr_data = data{idx}.Sg;
    hold on
    for j = 1:Npools
        fn = strrep(do_cells{j},'-',''); %structure fieldname for pool data
        plot(curr_data.(fn),'LineWidth',lnsz)
    end
end
Nobs = numel(curr_data.(fn));
xlim([1,Nobs])
ylabel(Ylab,'fontweight','b')
xlabel('seconds','FontWeight','b')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
legend(leg_labs,'Box','off')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'synaptic_gating'),'-djpeg')





%2) firing rate 
Ylab = {'E-stay','spike rate (Hz)'};
figure
for idx = 1:Njobs
    curr_data = data{idx}.Sr;
    hold on
    
    for j = 1:Npools
        fn = strrep(do_cells{j},'-',''); %structure fieldname for pool data
        plot(curr_data.(fn),'LineWidth',lnsz)
    end
end
Nobs = numel(curr_data.(fn));
xlim([1,Nobs])
ylabel(Ylab,'fontweight','b')
xlabel('seconds','FontWeight','b')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
legend(leg_labs,'Box','off')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'spikerates'),'-djpeg')



%3) depression
Ylab = {'E-stay','synaptic depression'};
figure
for idx = 1:Njobs
    curr_data = data{idx}.Dmu_fast;
    hold on
    for j = 1:Npools
        fn = strrep(do_cells{j},'-',''); %structure fieldname for pool data
        plot(curr_data.(fn),'LineWidth',lnsz)
    end
end
Nobs = numel(curr_data.(fn));
xlim([1,Nobs])
ylabel(Ylab,'fontweight','b')
xlabel('seconds','FontWeight','b')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
legend(leg_labs,'Box','off')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'depression'),'-djpeg')

%4) rate x depression
Ylab = {'E-stay','spikerate x depression'};
figure
for idx = 1:Njobs
    Sr = data{idx}.Sr;
    D = data{idx}.Dmu_fast;
    hold on
    
    for j = 1:Npools
        fn = strrep(do_cells{j},'-',''); %structure fieldname for pool data
        %rate time depression
        curr_data = Sr.(fn) .* D.(fn); 
        plot(curr_data,'LineWidth',lnsz)
    end
end
Nobs = numel(curr_data);
xlim([1,Nobs])
ylabel(Ylab,'fontweight','b')
xlabel('seconds','FontWeight','b')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
legend(leg_labs,'Box','off')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'Sr_X_depression'),'-djpeg')









function data = get_data(options)

%get this stuff
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%just look at some stuff here
savename = fullfile(options.save_dir,options.sim_name);
load(savename,'options','sim_results')
load(sprintf('%s_D',savename),'Drec_fast','Drec_slow')
load(sprintf('%s_V',savename),'Vrec')
load(sprintf('%s_spikes',savename),'spikes')
load(sprintf('%s_S',savename),'Srec')

term_idx = isnan(Vrec(1,:)); %when simulation terminated (beginning of 3rd stay state)
term_idx = find(term_idx,1,'first');
if ~isempty(term_idx)
    %now truncate
    spikes = spikes(:,1:term_idx);
    Drec_fast = Drec_fast(:,1:term_idx);
    Drec_slow = Drec_slow(:,1:term_idx);
    Vrec = Vrec(:,1:term_idx);
    Srec = Srec(:,1:term_idx);
end

timestep = options.timestep;


window_sz = 100e-3; %for moving average

data.Sr = sim_windowrate(spikes,timestep,celltype,window_sz);

data.Dmu_fast = simple_pool_avg(Drec_fast,celltype);
data.Dmu_slow = simple_pool_avg(Drec_slow,celltype);

data.Sg = sim_windowrate(Srec,timestep,celltype,window_sz);
data.Sg = structfun(@(x) x.* timestep ,data.Sg,'UniformOutput',false); %undo hz conversion
data.options = options;

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
end