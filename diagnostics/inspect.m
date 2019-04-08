clear
clc
format compact
hold off;close all

%investigating model behavior

addpath('../')
jobID = 4;
sname = 'diag_newswitch_crit';
%jobID = str2num(getenv('JID'));
%sname = getenv('SIM_NAME'); %'diag_EtoIfixed'
%my model
%---setup---------------------
options = set_options('modeltype','diagnostics','comp_location','woodstock','sim_name',sname,'jobID',jobID);


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
%lol = sim_results{4}

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


lnsz = 3; %spikerate plots
fontsz = 12;


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


Dmu_fast = simple_pool_avg(Drec_fast,celltype);
Dmu_slow = simple_pool_avg(Drec_slow,celltype);
valid_Drange = @(x) all(structfun(@max,x) < 1+eps) && all(structfun(@max,x) > 0-eps);

figure;
%plot the aggregated timecourses
ax(1) = subplot(2,1,1);hold on
plot(Dmu_fast.Estay,'Linewidth',lnsz)
plot(Dmu_fast.Eswitch,'Linewidth',lnsz)
plot(Dmu_fast.Istay,'Linewidth',lnsz)
plot(Dmu_fast.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Fast','depression'},'FontWeight','b')
%xlabel('time (s)')
set(gca,'FontSize',fontsz);axis tight;hold off
if valid_Drange(Dmu_fast),ylim([0,1]);end
axP = get(gca,'Position');
[lp, ~] = legend({'E-stay','E-switch','I-stay','I-switch'},'FontWeight','b',...
    'Location','northoutside','Box','off','Orientation','horizontal');
set(gca, 'Position', axP)
lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];
ax(2) = subplot(2,1,2);hold on
plot(Dmu_slow.Estay,'Linewidth',lnsz)
plot(Dmu_slow.Eswitch,'Linewidth',lnsz)
plot(Dmu_slow.Istay,'Linewidth',lnsz)
plot(Dmu_slow.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Slow','depression'},'FontWeight','b')
xlabel('time (s)')
set(gca,'FontSize',fontsz);axis tight;hold off
if valid_Drange(Dmu_slow),ylim([0,1]);end
linkaxes(ax,'xy')
print(fullfile(fig_dir,'depression'),'-djpeg')


figure;
%plot fast only depression
plot(Dmu_fast.Estay,'Linewidth',lnsz);hold on
plot(Dmu_fast.Eswitch,'Linewidth',lnsz)
plot(Dmu_fast.Istay,'Linewidth',lnsz)
plot(Dmu_fast.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Fast','depression'},'FontWeight','b')
xlabel('time (s)')
title('Depression')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
if valid_Drange(Dmu_fast),ylim([0,1]);end
print(fullfile(fig_dir,'depression_fastonly'),'-djpeg')


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
set(gca,'FontSize',fontsz);axis tight;hold off
print(fullfile(fig_dir,'spikerates'),'-djpeg')


figure;
%plot on/off pools seperately 
[on_range,off_range] = get_state_limits(S);
subplot(2,1,1);hold on
plot(S.Estay,'Linewidth',lnsz)
plot(S.Eswitch,'Linewidth',lnsz)
plot(S.Istay,'Linewidth',lnsz)
plot(S.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Inactive pool Hz','(75ms bins)'},'FontWeight','b');
set(gca,'FontSize',fontsz);axis tight;hold off
ylim(on_range)
axP = get(gca,'Position');
[lp, ~] = legend({'E-stay','E-switch','I-stay','I-switch'},'FontWeight','b',...
    'Location','northoutside','Box','off','Orientation','horizontal');
set(gca, 'Position', axP)
lp.Position = [(1-lp.Position(3))/2,1-lp.Position(4),lp.Position(3:4)];
subplot(2,1,2);hold on
plot(S.Estay,'Linewidth',lnsz)
plot(S.Eswitch,'Linewidth',lnsz)
plot(S.Istay,'Linewidth',lnsz)
plot(S.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
ylabel({'Inactive pool Hz','(75ms bins)'},'FontWeight','b');xlabel('time (s)')
set(gca,'FontSize',fontsz);axis tight;hold off
ylim(off_range)
print(fullfile(fig_dir,'spikerates_split'),'-djpeg')


figure;
plot(S.Estay-S.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({'Mean pool Hz  (75ms bins)'})
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
print(fullfile(fig_dir,'spikerates_Edifference'),'-djpeg')


Sg = simple_pool_avg(Srec,celltype);

%plot the aggregated timecourses
figure;
hold on
plot(Sg.Estay,'Linewidth',lnsz)
plot(Sg.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel('Mean pool S-gating')
xlabel('time (s)')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
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
ylabel('Mean pool S-gating')
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off %;ylim([0,.17])
print(fullfile(fig_dir,'syn_gating_allcells'),'-djpeg')
savefig(gcf(),fullfile(fig_dir,'syn_gating_allcells'),'compact')


%plot the aggregated timecourses
figure;
plot(Sg.Estay-Sg.Eswitch,'Linewidth',lnsz);hold on
plot(Sg.Istay-Sg.Iswitch,'Linewidth',lnsz);hold off
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel('Mean pool S-gating')
xlabel('time (s)')
legend({'E-stay minus E-switch','I-stay minus I-switch'},'location','northoutside','Orientation','horizontal')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
print(fullfile(fig_dir,'syn_gating_allcells_difference'),'-djpeg')

figure
%plot(Sg.Estay-Sg.Eswitch,'Linewidth',lnsz)
%Xticks = num2cell(get(gca,'Xtick'));
%Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
%set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
plot(0:timestep:options.tmax, Sg.Estay-Sg.Eswitch,'Linewidth',lnsz)
above_state_thresh = Sg.Estay-Sg.Eswitch;
above_state_thresh(abs(above_state_thresh) < options.state_test_thresh) = NaN;
hold on;plot(0:timestep:options.tmax, above_state_thresh,'Linewidth',lnsz);hold off
title('Synaptic gating')
ylabel('Mean pool S-gating')
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
print(fullfile(fig_dir,'syn_gating_difference'),'-djpeg')


%I want to see how long it takes to register a switch without "undecided" states
above_state_thresh = abs(Sg.Estay-Sg.Eswitch) > options.state_test_thresh;
above_state_thresh = diff(above_state_thresh) == -1;
tvec = timestep:timestep:options.tmax; %start with 1 b/c of diff
above_state_thresh = tvec(above_state_thresh);
above_state_thresh = diff(above_state_thresh);
above_state_thresh = above_state_thresh(above_state_thresh < .4); %doesn't take a half second
above_state_thresh = above_state_thresh(above_state_thresh > 25e-3); %takes longer than 25 ms
figure;histogram(above_state_thresh ./ 1e-3,numel(above_state_thresh));xlabel('transition detection time ms')
hold on 
above_state_thresh = Sg.Eswitch-Sg.Estay > options.state_test_thresh;
above_state_thresh = diff(above_state_thresh) == -1;
tvec = timestep:timestep:options.tmax; %start with 1 b/c of diff
above_state_thresh = tvec(above_state_thresh);
above_state_thresh = diff(above_state_thresh);
above_state_thresh = above_state_thresh(above_state_thresh < .3); %doesn't take a half second
above_state_thresh = above_state_thresh(above_state_thresh > 25e-3); %takes longer than 25 ms


window_sz = 100e-3;
Sg_smooth = sim_windowrate(Srec,timestep,celltype,window_sz);
Sg_smooth = structfun(@(x) x.* timestep ,Sg_smooth,'UniformOutput',false); %undo hz conversion

figure;
hold on
plot(Sg_smooth.Estay,'Linewidth',lnsz)
plot(Sg_smooth.Eswitch,'Linewidth',lnsz)
plot(Sg_smooth.Istay,'Linewidth',lnsz)
plot(Sg_smooth.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Smoothed synaptic gating')
ylabel('smoothed S-gating')
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off ;ylim([0,.17])
print(fullfile(fig_dir,'syn_gating_allcells_smoothed'),'-djpeg')


durations = sim_results{1};
if numel(durations) > 1
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
    timecourse = [timecourse(:,1),num2cell(cellfun(@(x,y) x-y,timecourse(:,1),timecourse(:,2))),timecourse(:,2:end)];
    %timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
    timecourse = cell2table(timecourse,'VariableNames',{'event_time','start','duration','sample_time','sample_number','state'});
    figure
    %state_durs = timecourse.duration(startsWith(timecourse.state,'stim'));
    state_durs = timecourse.duration(~strcmpi(timecourse.state,'undecided')); %leave states included 
    if numel(state_durs) > 0
        if numel(state_durs) < 15
            histogram(state_durs,numel(state_durs))
        else
            histogram(state_durs)
        end
        %title(sprintf('stay-state durations\nmu = %.2f',mean(state_durs)))
        %ylabel('seconds');xlabel('frequecy');set(gca,'FontSize',fontsz)
        title(sprintf('Network:    E-I = %.2f,    I-E = %.2f\nmu duration = %.2fs',options.EtoI,options.ItoE,mean(state_durs)))
        ylabel('state duration (s)');xlabel('frequecy');set(gca,'FontSize',fontsz)
        print(fullfile(fig_dir,'state_durations'),'-djpeg')
    end
    %add paramters, print
end

figure
title(sprintf('Network:    E-I = %.2f,    I-E = %.2f',options.EtoI,options.ItoE))
set(gca,'FontSize',fontsz)
print(fullfile(fig_dir,'parameter_title'),'-djpeg')


% %membrane voltage distributions
% Ve = Vrec(celltype.excit,:);
% Ve = Ve(:) ./ 1e-3; %convert to mV units 
% Vi = Vrec(celltype.inhib,:);
% Vi = Vi(:) ./ 1e-3; %convert to mV units 
% 
% figure;orient tall
% subplot(2,1,1)
% histogram(Ve)
% title(sprintf('Excitatory cells (dt = %.2f ms)',timestep / 1e-3))
% legend(sprintf('\\mu = %.1f mV\\newline\\sigma = %1.f mV',mean(Ve),std(Ve)))
% set(gca,'FontSize',fontsz,'FontWeight','b')
% subplot(2,1,2)
% histogram(Vi)
% title(sprintf('Inhibitory cells (dt = %.2f ms)',timestep / 1e-3))
% legend(sprintf('\\mu = %.1f mV\\newline\\sigma = %1.f mV',mean(Vi),std(Vi)))
% set(gca,'FontSize',fontsz,'FontWeight','b')
% xlabel('membrane potential (mV)')
% print(fullfile(fig_dir,'Vm_dists'),'-djpeg')



%for printing 
TCfile = cellfun(@(x) sprintf('%.3f',x),table2cell(timecourse(:,1:end-1)),'UniformOutput',false); %for printing
TCfile = [TCfile,timecourse.state];
TCfile = cell2table(TCfile,'VariableNames',{'event_time','start','duration','sample_time','sample_number','state'});
writetable(TCfile,fullfile(fig_dir,'event_info.txt'),'Delimiter','|')

% %this is what happened, and when
% %{timeidx,statecount,sample_clock,stim_label}
% durations = sim_results{1};
% timecourse = size(durations);
% timecourse(2) = timecourse(2) + 1; 
% timecourse = cell(timecourse);
% timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
% timecourse(:,3) = cellfun(@(x) mod(x,sum(options.stim_pulse)),timecourse(:,3),'UniformOutput',false);
% %current sample's onset, rounding is needed for subsequent operations 
% timecourse(:,4) = cellfun(@(x,y) round(x-y,2),timecourse(:,1),timecourse(:,3),'UniformOutput',false);
% samp_onsets = unique(cat(1,timecourse{:,4})); %like unique won't work properly here without rounding 
% timecourse(:,4) = cellfun(@(x) find(x==samp_onsets),timecourse(:,4),'UniformOutput',false); %would also break without rounding
% timecourse(:,end) = durations(:,end);
% timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
% timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});
% writetable(timecourse,fullfile(fig_dir,'event_info.txt'),'Delimiter','|')
% 
% %figure out the difference between end of stimulus & next stimulus onset
% [out,info] = find_stay_durations(durations,options,'verify');
% onset_diff = NaN(size(info,1),1);
% time_inds = cell2mat(cellfun(@(x) x*options.timestep,durations(:,1),'UniformOutput',false));
% u2_times = NaN(size(onset_diff)); %time for second undecided stte 
% for idx = 1:size(info,1)
%     Tend = info{idx,1};
%     curr_stim = find(time_inds == Tend);
%     if curr_stim+4 <= size(timecourse,1)
%         curr_stim = timecourse(curr_stim:curr_stim+4,:);
%         if ~strcmpi(curr_stim{3,end},'leave') || ~startsWith(curr_stim{end,end},'stim')
%             fprintf('problem with item #%i\r',idx) %almost certainly a "stim-leave-leave" sequence, no biggie
%         else
%             onset_diff(idx) = str2num(curr_stim.event_time{end-1}) - str2num(curr_stim.event_time{1});
%             u2_times(idx) = str2num(curr_stim.duration{4});
%         end
%     end
% end
% Tinvalid = isnan(onset_diff);
% fprintf('valid transitions = %i/%i\n',sum(~Tinvalid),numel(Tinvalid))
% if nansum(onset_diff) > 0
%     onset_diff = onset_diff(~Tinvalid);
%     u2_times = u2_times(~Tinvalid);
%     figure()
%     histogram(onset_diff,numel(onset_diff));hold on
%     histogram(u2_times,numel(u2_times));hold off
%     title('stim onset differences')
%     legend('total transition','undecided #2','Location','best')
% end



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

%spikerates by sliding window...

% window_sz = 50e-3;
% S = sim_windowrate(spikes,timestep,celltype,window_sz);
% 
% figure;
% %plot the aggregated timecourses
% hold on
% plot(S.Estay,'Linewidth',lnsz)
% plot(S.Eswitch,'Linewidth',lnsz)
% %plot(S.Istay,'Linewidth',lnsz)
% %plot(S.Iswitch,'Linewidth',lnsz)
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% title('Spiking')
% ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
% xlabel('time (s)')
% %legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
% legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
% set(gca,'FontSize',fontsz);axis tight;hold off
% print(fullfile(fig_dir,'spikerates_slidingwin'),'-djpeg')
% savefig(gcf(),fullfile(fig_dir,'spikerates_slidingwin'),'compact')
% 
% figure;
% %plot the aggregated timecourses
% hold on
% plot(S.Estay,'Linewidth',lnsz)
% plot(S.Eswitch,'Linewidth',lnsz)
% plot(S.Istay,'Linewidth',lnsz)
% plot(S.Iswitch,'Linewidth',lnsz)
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% title('Spiking')
% ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
% xlabel('time (s)')
% legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
% set(gca,'FontSize',fontsz);axis tight;hold off
% print(fullfile(fig_dir,'spikerates_slidingwin_allcells'),'-djpeg')
% savefig(gcf(),fullfile(fig_dir,'spikerates_slidingwin_allcells'),'compact')
