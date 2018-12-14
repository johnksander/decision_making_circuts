clear
clc
format compact
hold off;close all

%investigating model behavior

addpath('../')

%my model
%---setup---------------------
t = 200;
options = set_options('modeltype','PS_stim',...
    'sim_name','timer',...
    'jobID',1,'tmax',t,'stim_pulse',[t,0],'stim_schedule','flexible',...
    'comp_location','woodstock','cut_leave_state',5e-3);


fig_dir = fullfile(options.save_dir,options.sim_name);
if ~isdir(fig_dir),mkdir(fig_dir);end

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
load(sprintf('%s_S',savename))

term_idx = isnan(Vrec(1,:)); %when simulation terminated (beginning of 3rd stay state)
term_idx = find(term_idx,1,'first');
if ~isempty(term_idx)
    %now truncate
    spikes = spikes(:,1:term_idx);
    Drec = Drec(:,1:term_idx);
    Vrec = Vrec(:,1:term_idx);
    Srec = Srec(:,1:term_idx);
end

timestep = .25e-3; %.25 milisecond timestep
lnsz = 3; %spikerate plots
fontsz = 12;

%this is what happened, and when
%{timeidx,statecount,sample_clock,stim_label}
durations = sim_results{1};
timecourse = size(durations);
timecourse(2) = timecourse(2) + 1; 
timecourse = cell(timecourse);
timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
timecourse(:,3) = cellfun(@(x) mod(x,sum(options.stim_pulse)),timecourse(:,3),'UniformOutput',false);
%current sample's onset, rounding is needed for subsequent operations 
timecourse(:,4) = cellfun(@(x,y) round(x-y,2),timecourse(:,1),timecourse(:,3),'UniformOutput',false);
samp_onsets = unique(cat(1,timecourse{:,4})); %like unique won't work properly here without rounding 
timecourse(:,4) = cellfun(@(x) find(x==samp_onsets),timecourse(:,4),'UniformOutput',false); %would also break without rounding
timecourse(:,end) = durations(:,end);
timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});
writetable(timecourse,fullfile(fig_dir,'event_info.txt'),'Delimiter','|')

%figure out the difference between end of stimulus & next stimulus onset
[out,info] = find_stay_durations(durations,options,'verify');
onset_diff = NaN(size(info,1),1);
time_inds = cell2mat(cellfun(@(x) x*options.timestep,durations(:,1),'UniformOutput',false));
u2_times = NaN(size(onset_diff)); %time for second undecided stte 
for idx = 1:size(info,1)
    Tend = info{idx,1};
    curr_stim = find(time_inds == Tend);
    if curr_stim+4 <= size(timecourse,1)
        curr_stim = timecourse(curr_stim:curr_stim+4,:);
        if ~strcmpi(curr_stim{3,end},'leave') || ~startsWith(curr_stim{end,end},'stim')
            fprintf('problem with item #%i\r',idx) %almost certainly a "stim-leave-leave" sequence, no biggie
        else
            onset_diff(idx) = str2num(curr_stim.event_time{end-1}) - str2num(curr_stim.event_time{1});
            u2_times(idx) = str2num(curr_stim.duration{4});
        end
    end
end
Tinvalid = isnan(onset_diff);
fprintf('valid transitions = %i/%i\n',sum(~Tinvalid),numel(Tinvalid))
if sum(onset_diff) > 0
    onset_diff = onset_diff(~Tinvalid);
    u2_times = u2_times(~Tinvalid);
    figure()
    histogram(onset_diff,numel(onset_diff));hold on
    histogram(u2_times,numel(u2_times));hold off
    title('stim onset differences')
    legend('total transition','undecided #2','Location','best')
end

figure()
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
ylabel({'Pool average  (75ms bins)'})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'depression'),'-djpeg')

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
ylabel({'Mean pool Hz  (75ms bins)'})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'spikerates'),'-djpeg')


%spikerates by sliding window...

window_sz = 50e-3;
S = sim_windowrate(spikes,timestep,celltype,window_sz);

figure;
%plot the aggregated timecourses
hold on
plot(S.Estay,'Linewidth',lnsz)
plot(S.Eswitch,'Linewidth',lnsz)
%plot(S.Istay,'Linewidth',lnsz)
%plot(S.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'spikerates_slidingwin'),'-djpeg')
savefig(gcf(),fullfile(fig_dir,'spikerates_slidingwin'),'compact')

figure;
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
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'spikerates_slidingwin_allcells'),'-djpeg')
savefig(gcf(),fullfile(fig_dir,'spikerates_slidingwin_allcells'),'compact')



figure;
plot(S.Estay-S.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'spikerates_slidingwin_difference'),'-djpeg')



%Sg = sim_windowrate(Srec,timestep,celltype,window_sz);
%Sg = structfun(@(x) x.* timestep ,Sg,'UniformOutput',false); %undo hz conversion
Sg.Estay = mean(Srec(celltype.excit & celltype.pool_stay,:),1);
Sg.Eswitch = mean(Srec(celltype.excit & celltype.pool_switch,:),1);
Sg.Istay = mean(Srec(celltype.inhib & celltype.pool_stay,:),1);
Sg.Iswitch = mean(Srec(celltype.inhib & celltype.pool_switch,:),1);

%plot the aggregated timecourses
figure;
hold on
plot(Sg.Estay,'Linewidth',lnsz)
plot(Sg.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'syn_gating'),'-djpeg')
savefig(gcf(),fullfile(fig_dir,'syn_gating'),'compact')

%plot the aggregated timecourses
figure;
hold on
plot(Sg.Estay,'Linewidth',lnsz)
plot(Sg.Eswitch,'Linewidth',lnsz)
plot(Sg.Istay,'Linewidth',lnsz)
plot(Sg.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
hold off
print(fullfile(fig_dir,'syn_gating_allcells'),'-djpeg')
savefig(gcf(),fullfile(fig_dir,'syn_gating_allcells'),'compact')



%plot the aggregated timecourses
figure;
plot(Sg.Estay-Sg.Eswitch,'Linewidth',lnsz);hold on
plot(Sg.Istay-Sg.Iswitch,'Linewidth',lnsz);hold off
% plot(Sg.Estay,'Linewidth',lnsz);hold on
% plot(Sg.Eswitch,'Linewidth',lnsz)
% plot(Sg.Istay,'Linewidth',lnsz)
% plot(Sg.Iswitch,'Linewidth',lnsz);hold off
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch','I-stay minus I-switch'},'location','northoutside','Orientation','horizontal')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)




figure
plot(Sg.Estay-Sg.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'syn_gating_difference'),'-djpeg')




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

function cell_data = sim_windowrate(cell_raster,timestep,celltype,window_sz)

num_binsamps = window_sz/timestep; %num samples in Xms
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




