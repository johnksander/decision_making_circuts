clear
clc
format compact
hold off;close all

%investigating model behavior

addpath('../')

%my model
%---setup---------------------
config_options.modeltype = 'PS_stim';
config_options.sim_name = 'diagnostics';
%specify network #1 slow w/ jobID
config_options.jobID = 2; %str2num(getenv('SGE_TASK_ID'));
config_options.tmax = 30; %trial simulation time (s)
config_options.force_back2stay = true;
config_options.stim_pulse = [1, 1]; %on, off (s)
options = set_options(config_options);

% %NOTE: I AM NOT removing the first artificial stay state here. 
% %---run-----------------------
% exit_status = false;
% while ~exit_status
%     [modelfile,exit_status] = diag_model(options);
% end
% %---cleanup-------------------
% driverfile = mfilename;
% backup_jobcode(options,driverfile,modelfile)
% delete(options.output_log) %no need for these right now


%get this stuff 
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
celltype = celltype_logicals(pool_options);


%just look at some stuff here 
savename = fullfile(options.save_dir,options.sim_name);
load(savename)
load(sprintf('%s_D',savename))
load(sprintf('%s_V',savename))
load(sprintf('%s_spikes',savename))

term_idx = isnan(Vrec(1,:)); %when simulation terminated (beginning of 3rd stay state)
term_idx = find(term_idx,1,'first');

%now truncate 
spikes = spikes(:,1:term_idx);
Drec = Drec(:,1:term_idx);
Vrec = Vrec(:,1:term_idx);

timestep = .25e-3; %.25 milisecond timestep
lnsz = 3; %spikerate plots
fontsz = 14;

%this is what happened, and when 
durations = sim_results{1};
timecourse = cell(size(cat(1,durations{:})));
timecourse(1:2:end,:) = durations{1};
timecourse(2:2:end,:) = durations{2};
timecourse(:,1) = cellfun(@(x) x*timestep,timecourse(:,1),'UniformOutput',false);
num_samps = cellfun(@(x) x./ sum(options.stim_pulse),timecourse(:,1),'UniformOutput',false);
event_times = cumsum(cell2mat(timecourse(:,1)));
timecourse = [num2cell(event_times),num_samps,timecourse];
timecourse(~cellfun(@isnan,timecourse(:,end)),end) = {'stay'};
timecourse(~cellfun(@ischar,timecourse(:,end)),end) = {'switch'};
timecourse = cell2table(timecourse,'VariableNames',{'event_time','samples','duration','state'});


%do the raster plot 
spikeplot = make_spikeplot(spikes);
%reorganize the cell groups 
raster = [spikeplot(celltype.pool_stay & celltype.excit,:);...
    spikeplot(celltype.pool_switch & celltype.inhib,:);...
    spikeplot(celltype.pool_switch & celltype.excit,:);...
    spikeplot(celltype.pool_stay & celltype.inhib,:)];
imagesc(raster)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
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
xlabel('time (s)','FontWeight','b')
title({'spikes','(spikes in matrix enlarged for visualization)'})
set(gca,'FontSize',fontsz)

print(fullfile(fig_dir,'raster'),'-djpeg')

%plot the aggregated timecourses
figure;

D = sim_spikerate(Drec,timestep,celltype);
D = structfun(@(x) x.* timestep ,D,'UniformOutput',false); %undo hz conversion

%plot the aggregated timecourses
hold on
plot(D.Estay,'Linewidth',lnsz)
plot(D.Eswitch,'Linewidth',lnsz)
plot(D.Istay,'Linewidth',lnsz)
plot(D.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Depression')
ylabel({'Pool average', '(75ms time-bins)'})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off

%spikerates
figure;

S = sim_spikerate(spikes,timestep,celltype);

%plot the aggregated timecourses
hold on
plot(S.Estay,'Linewidth',lnsz)
plot(S.Eswitch,'Linewidth',lnsz)
plot(S.Istay,'Linewidth',lnsz)
plot(S.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({'Mean pool Hz','(75ms time-bins)'})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off





function cell_data = sim_spikerate(cell_raster,timestep,celltype)

num_binsamps = 75e-3/timestep; %num samples in 2ms
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
