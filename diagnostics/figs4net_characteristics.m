clear
clc
format compact
hold off;close all

%show network characteristics in four paneled plot
%this requires saved datafiles from a "diagnostic" model run


addpath('../')

%specify what results
Sname = 'example_behavior';
jobID = 444;
options = set_options('modeltype','diagnostics','comp_location','bender','sim_name',Sname,'jobID',jobID);

%figure options
lnsz = 3;
fontsz = 12;
shift_start = 2; %plot starting x seconds into simulation data 

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

if shift_start > 0
    start_ind = round(shift_start / timestep);
    %now truncate
    spikes = spikes(:,start_ind:end);
    Drec_fast = Drec_fast(:,start_ind:end);
    Drec_slow = Drec_slow(:,start_ind:end);
    Vrec = Vrec(:,start_ind:end);
    Srec = Srec(:,start_ind:end);
end


%we want a raster, D-fast, rate plot, D-slow (in subplot order)

figure;subplot(2,2,1) %do the raster plot
 
durr_spike = 20e-3; %how wide is each spike (duration, in seconds ) 
raster = get_raster_img(spikes,timestep,durr_spike);

%reorganize the cell groups
raster = [raster(celltype.pool_stay & celltype.excit,:);...
    raster(celltype.pool_switch & celltype.inhib,:);...
    raster(celltype.pool_switch & celltype.excit,:);...
    raster(celltype.pool_stay & celltype.inhib,:)];
imagesc(raster)

colormap(flipud(colormap('gray')));

Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);

Yticks = num2cell(get(gca,'YLim'));
Yticks = Yticks{2};
Yticks = [Yticks*.25,Yticks*.75];
Ylabs = {'stay pool','leave pool'};
ytickangle(90)
set(gca,'Xdir','normal','Ytick',Yticks,'YTickLabel', Ylabs);
% set(gca,'Xdir','normal','Ytick',Yticks,'YTickLabel', {''});
% text(-2e3,Yticks*.2,'E-stay','HorizontalAlignment','center','Rotation',90)
% text(-2e3,Yticks*.45,'I-stay','HorizontalAlignment','center','Rotation',90)
% text(-2e3,Yticks*.7,'E-switch','HorizontalAlignment','center','Rotation',90)
% text(-2e3,Yticks*.95,'I-switch','HorizontalAlignment','center','Rotation',90)
xlabel('seconds')
%title({'spikes','(spikes in matrix enlarged for visualization)'})
set(gca,'FontSize',fontsz)

%print(fullfile(fig_dir,'raster'),'-djpeg')



%spikerates
subplot(2,2,3)

window_sz = 100e-3;
S = sim_windowrate(spikes,timestep,celltype,window_sz);

%plot the aggregated timecourses
hold on
plot(S.Estay,'Linewidth',lnsz)
plot(S.Eswitch,'Linewidth',lnsz)
plot(S.Istay,'Linewidth',lnsz)
plot(S.Iswitch,'Linewidth',lnsz);%axis tight
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Mean pool spiking (Hz)'})
xlabel('seconds')
set(gca,'FontSize',fontsz)
axP = get(gca,'Position');
[lp, ~] = legend({'E-stay','E-switch','I-stay','I-switch'},'FontWeight','b',...
    'Location','northoutside','Box','off','Orientation','horizontal','FontSize',fontsz);
set(gca, 'Position', axP)
lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];

%one big legend over the whole joint.. 





Dmu_fast = simple_pool_avg(Drec_fast,celltype);
Dmu_slow = simple_pool_avg(Drec_slow,celltype);
valid_Drange = @(x) all(structfun(@max,x) < 1+eps) && all(structfun(@max,x) > 0-eps);

%plot the aggregated timecourses
ax(1) = subplot(2,2,2);hold on
plot(Dmu_fast.Estay,'Linewidth',lnsz)
plot(Dmu_fast.Eswitch,'Linewidth',lnsz)
plot(Dmu_fast.Istay,'Linewidth',lnsz)
plot(Dmu_fast.Iswitch,'Linewidth',lnsz);
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
Xlabs = strrep(Xlabs,'.0','');
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Fast','depression'},'FontWeight','b')
xlabel('seconds')
set(gca,'FontSize',fontsz);hold off; %axis tight
%if valid_Drange(Dmu_fast),ylim([0,1]);end
ax(2) = subplot(2,2,4);hold on
plot(Dmu_slow.Estay,'Linewidth',lnsz)
plot(Dmu_slow.Eswitch,'Linewidth',lnsz)
plot(Dmu_slow.Istay,'Linewidth',lnsz)
plot(Dmu_slow.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
Xlabs = strrep(Xlabs,'.0','');
ylabel({'Slow','depression'},'FontWeight','b')
xlabel('seconds')
set(gca,'FontSize',fontsz); hold off; %axis tight
%if valid_Drange(Dmu_slow),ylim([0,1]);end
linkaxes(ax,'xy')


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

