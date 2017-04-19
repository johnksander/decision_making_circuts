function behavior_test_JK_noD(current_weight)

% clear
% clc
close all
% format compact

%original code: switching_simulationv16_Eonlyv2_noD.m

basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
addpath(fullfile(basedir,'helper_functions'))
testdir = fullfile(basedir,'inhibition_test');
figdir =  fullfile(testdir,'sim_figures');
if ~isdir(figdir)
    mkdir(figdir)
end

%current_weight = 1; %basecurr modifier
fileID = strrep(num2str(current_weight),'.','p');
simname = ['JK_' fileID];
savedir = fullfile(figdir,simname);
if ~isdir(savedir)
    mkdir(savedir)
end

output_vars = 'on';
switch output_vars
    case 'on'
        vardir = fullfile(testdir,'output_vars');
        if ~isdir(vardir)
            mkdir(vardir)
        end
end


%same as v16_Eonlyv2, but I'm removing synaptic depression here.

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 200;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.5 .5]; %proportion excitable % inhibitory
pool_options.p_conn = .5; %connection probability 50%
%--------------------------------------------------------------------------
%build the connectivity matrix
%make celltype logicals
celltype = celltype_logicals(pool_options);
%make these vectors into logical matricies
[EtoE,EtoI,ItoE] = connection_logicals(celltype,pool_options.num_cells);
%make connection scheme based off connection matricies
connection_scheme = EtoE | EtoI | ItoE;  %plan is:  EtoE | EtoI | ItoE;
%make synaptic connection matrix
W = rand(pool_options.num_cells) < pool_options.p_conn; %connection matrix
W = double(W & connection_scheme); %filter out connections not allowed by scheme
%modify connection weights
W(W > 0 & EtoE) = .25 /12;
W(W > 0 & ItoE) = 5 /7;
W(W > 0 & EtoI) = .075 /7; %changed since v16
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
current_val = .1e-9; %applied current in nano amps %current & noise adjusted from v16
current_val = current_val * current_weight;
noise_sigma = 10e-12; %noise variance in picoamps
current_pulse = 'off'; %switch for adding transient pulses
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential, excitatory
Erev(celltype.inhib) = -70e-3; %reversal potential, inhibitory
Gg = NaN(pool_options.num_cells,1); %max connductance
Gg(celltype.excit) = 10e-9; %excitatory max connductance microSiemens
Gg(celltype.inhib) = 10e-9; %inhibitory max connductance microSiemens
Pr = NaN(pool_options.num_cells,1); %release probability
% Pr(celltype.excit) = .2; %excitatory release probability
% Pr(celltype.inhib) = .2; %inhibitory release probability
Pr(celltype.excit) = 1; %excitatory release probability
Pr(celltype.inhib) = 1; %inhibitory release probability
Td = 5e-3;%synaptic depression time constant, seconds
%Td = .3;
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50e-3; %excitatory gating time constant, ms
Tsyn(celltype.inhib) = 10e-3; %inhibitory gating time constant, ms
%----cell basics--------------
El = -70e-3; %leak potential mV
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%---adaptation current--------
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
a = 1.25e-9; %adaptation conductance, in nano Siemens
b = NaN(pool_options.num_cells,1); %increase adaptation current by max value, nano amps
b(celltype.excit) = .25e-9; %excitatory
b(celltype.inhib) = .25e-9; %inhibitory
Tsra = NaN(pool_options.num_cells,1);%adaptation current time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
%----timecourse---------------
tmax = 100; %simulation end (s)
timestep = .25e-3; %.25 milisecond timestep
timevec = 0:timestep:tmax;
switch current_pulse
    case 'on'
        ptime_stay = [.5,.55];
        ptime_switch = [1.5,1.55];
        
        pulse_current = .1e-9;
        pulse_params = {celltype.pool_stay,pulse_current,ptime_stay(1),ptime_stay(2);...
            celltype.pool_switch,pulse_current,ptime_switch(1),ptime_switch(2);...
            celltype.pool_stay,pulse_current,2.5,2.55;...
            celltype.pool_stay,pulse_current,8,8.05;...
            celltype.pool_switch,pulse_current,12.1,12.15}; %logical, current, start, end
        
        %pulse_params = {celltype.pool_stay,.001e-9,0,25};
end
%-------------------------------------------------------------------------
%preallocate variables

V = NaN(pool_options.num_cells,numel(timevec)); %membrane potential
V(:,1) = El; %inital value of membrane potential is leak potential

Iapp = zeros(size(V)); %applied current vector

%ONLY EXCITATORY
Iapp(celltype.excit,:) = Iapp(celltype.excit,:) + current_val; %ONLY EXCITATORY
%ONLY EXCITATORY

noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment
Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise
switch current_pulse
    case 'on'
        Iapp = addpulse(pulse_params,Iapp,timevec); %add pulses
end

Isra = NaN(size(V)); %adaptation current
Isra(:,1) = 0; %initialize at zero??

Sg = NaN(size(V)); %synaptic gating
Sg(:,1) = 0; %initalize at zero??
D = NaN(size(V)); %synaptic depression
D(:,1) = 1; %initalize at one

spikes = zeros(size(V));
durations = cell(1,2); %record duration times (stay, switch)
state.now = logical([1 0]);
%pick one state to start with, guess it doesn't really matter and
%avoids an if-statement check on every iteration (if initalized as NaN or something)
state.count = 0;
state.pools2compare = [celltype.pool_stay,celltype.pool_switch]; %pass in this format, avoid many computations

for idx = 2:numel(timevec)
    
    
    I = Iapp(:,idx-1);
    I = I + (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
    I = I + (unique(Erev(celltype.inhib)) - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*unique(Gg(celltype.inhib));
    dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm) - Isra(:,idx-1) + I;
    V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
    
    dt_Isra = ((a .* (V(:,idx-1) - El)) - Isra(:,idx-1)) ./ Tsra ; %adaptation current rate of change
    Isra(:,idx) = (dt_Isra .* timestep) + Isra(:,idx-1);
    %D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
    D(:,idx) = 1;
    Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
    
    if  sum(V(:,idx) > spike_thresh) > 0
        spiking_cells = V(:,idx) > spike_thresh;
        Isra(spiking_cells,idx) = Isra(spiking_cells,idx) + b(spiking_cells); %increase adaptation current by max value
        Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
            (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
        D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
        V(spiking_cells,idx) = Vreset;
        spikes(spiking_cells,idx) = 1;
    end
    D(:,idx) = 1;
    %[state,durations] = test4switch(Sg(:,idx),state,durations);
    
end

example_cells = [85 20];
% 
% figure(1)
% spikeplot = make_spikeplot(spikes);
% imagesc(spikeplot)
% Xtick_seconds=0:.5:tmax;
% Xtick_pos = find(mod(timevec,.5) == 0);
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% title({'spikes','(spikes in matrix enlarged for visualization)'})
% ylabel('cell')
% xlabel('time (s)')
% print(fullfile(savedir,['spikeplot_' fileID]),'-djpeg')
% 
% figure(2)
% imagesc(Iapp)
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% title('current')
% ylabel('cell')
% xlabel('time (s)')
% print(fullfile(savedir,['current_' fileID]),'-djpeg')
% 
% figure(3)
% subplot(2,1,1)
% plot(D(example_cells(1),:))
% hold on
% plot(D(example_cells(2),:))
% hold off
% title('synaptic depression (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% subplot(2,1,2)
% plot(D(example_cells(1) + 100,:))
% hold on
% plot(D(example_cells(2) + 100,:))
% hold off
% title('synaptic depression (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% print(fullfile(savedir,['synaptic_depression_' fileID]),'-djpeg')
% 
% 
% 
% 
% figure(4)
% subplot(2,1,1)
% plot(Sg(example_cells(1),:))
% hold on
% plot(Sg(example_cells(2),:))
% title('synaptic gating (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% subplot(2,1,2)
% plot(Sg(example_cells(1) + 100,:))
% hold on
% plot(Sg(example_cells(2) + 100,:))
% title('synaptic gating (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% print(fullfile(savedir,['synaptic_gating_' fileID]),'-djpeg')
% 
% 
% figure(5)
% subplot(2,1,1)
% plot(V(example_cells(1),:))
% hold on
% plot(V(example_cells(2),:))
% title('membrane voltage (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% subplot(2,1,2)
% plot(V(example_cells(1) + 100,:))
% hold on
% plot(V(example_cells(2) + 100,:))
% title('membrane voltage (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% print(fullfile(savedir,['Vm_' fileID]),'-djpeg')
% 
% 
% figure(6)
% subplot(2,1,1)
% plot(Isra(example_cells(1),:))
% hold on
% plot(Isra(example_cells(2),:))
% hold off
% legend('inhibitory','excitatory')
% title('Adaptation current (stay pool example cells)')
% xlabel('time (s)')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% subplot(2,1,2)
% plot(Isra(example_cells(1) + 100,:))
% hold on
% plot(Isra(example_cells(2) + 100,:))
% hold off
% legend('inhibitory','excitatory')
% title('Adaptation current (switch pool example cells)')
% xlabel('time (s)')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% print(fullfile(savedir,['adaptation_current_' fileID]),'-djpeg')
% 
% 
% figure(7)
% subplot(2,1,1)
% Irateplot = make_spikerate_plot(spikes,celltype.pool_stay & celltype.inhib,timestep);
% Irateplot(isnan(Irateplot)) = 0;
% plot(Irateplot,'linewidth',2)
% hold on
% Erateplot = make_spikerate_plot(spikes,celltype.pool_stay & celltype.excit,timestep);
% Erateplot(isnan(Erateplot)) = 0;
% plot(Erateplot,'linewidth',2)
% title('stay pool cells')
% legend('inhibitory','excitatory')
% ylabel({'average cell rate (Hz)','(mean cell 1/ISI at each timepoint)'})
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% xlabel('time (s)')
% hold off
% subplot(2,1,2)
% Irateplot = make_spikerate_plot(spikes,celltype.pool_switch & celltype.inhib,timestep);
% Irateplot(isnan(Irateplot)) = 0;
% plot(Irateplot,'linewidth',2)
% hold on
% Erateplot = make_spikerate_plot(spikes,celltype.pool_switch & celltype.excit,timestep);
% Erateplot(isnan(Erateplot)) = 0;
% plot(Erateplot,'linewidth',2)
% title('switch pool cells')
% legend('inhibitory','excitatory')
% ylabel({'average cell rate (Hz)','(mean cell 1/ISI at each timepoint)'})
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% xlabel('time (s)')
% hold off
% print(fullfile(savedir,['spikerates_' fileID]),'-djpeg')
% 
% 
% figure(8)
% cell2plot = find(celltype.excit,1,'first');
% histogram(Iapp(cell2plot,:))
% title(['current distribution, cell #' num2str(cell2plot)])
% print(fullfile(savedir,['current_dist_excit_' fileID]),'-djpeg')
% 
% figure(9)
% cell2plot = find(celltype.inhib,1,'first');
% histogram(Iapp(cell2plot,:))
% title(['current distribution, cell #' num2str(cell2plot)])
% print(fullfile(savedir,['current_dist_inhib_' fileID]),'-djpeg')
% 
% 
% 
% % disp('---during pulse:')
% % pulse_time = timevec >= ptime(1) & timevec <= ptime(2);
% % Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
% % Erate = mean(Erate) / (sum(pulse_time) * timestep);
% % disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
% % Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
% % Irate = mean(Irate) / (sum(pulse_time) * timestep);
% % disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
% disp('---during t+ two seconds (stay):')
% pulse_time = timevec >= 2;
% Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
% Erate = mean(Erate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
% Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
% Irate = mean(Irate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
% disp('---during t+ two seconds (switch):')
% pulse_time = timevec >= 2;
% Erate = sum(spikes(celltype.pool_switch & celltype.excit,pulse_time),2);
% Erate = mean(Erate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
% Irate = sum(spikes(celltype.pool_switch & celltype.inhib,pulse_time),2);
% Irate = mean(Irate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))



switch output_vars
    case 'on'
        muEstay_Sg = mean(Sg(celltype.pool_stay & celltype.excit,:)); %Sg means for both both pools
        muEswitch_Sg = mean(Sg(celltype.pool_switch & celltype.excit,:));
        
        Erate_stay = make_spikerate_plot(spikes,celltype.pool_stay & celltype.excit,timestep);
        Erate_stay(isnan(Erate_stay)) = 0; %stay
        Erate_switch = make_spikerate_plot(spikes,celltype.pool_switch & celltype.excit,timestep);
        Erate_switch(isnan(Erate_switch)) = 0; %switch
        
        Irate_stay = make_spikerate_plot(spikes,celltype.pool_stay & celltype.inhib,timestep);
        Irate_stay(isnan(Irate_stay)) = 0; %stay
        Irate_switch = make_spikerate_plot(spikes,celltype.pool_switch & celltype.inhib,timestep);
        Irate_switch(isnan(Irate_switch)) = 0; %switch...  angry variables, they're irate!! 
        
        save(fullfile(vardir,simname),'muEstay_Sg','muEswitch_Sg',...
            'Erate_stay','Erate_switch','Irate_stay','Irate_switch')
end




% figure(2)
% histogram(durations{1} * timestep)
% hold on
% histogram(durations{2} * timestep)
% legend('stay','switch')
% title('state durations')
% ylabel('frequency')
% xlabel('time (s)')
% hold off

%mean plots:

%
% figure(4)
% subplot(2,1,1)
% plot(mean(Sg(celltype.pool_stay & celltype.inhib,:)))
% hold on
% plot(mean(Sg(celltype.pool_stay & celltype.excit,:)))
% title('synaptic gating (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% subplot(2,1,2)
% plot(mean(Sg(celltype.pool_switch & celltype.inhib,:)))
% hold on
% plot(mean(Sg(celltype.pool_switch & celltype.excit,:)))
% title('synaptic gating (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
%
%
% figure(5)
% subplot(2,1,1)
% plot(mean(V(celltype.pool_stay & celltype.inhib,:)))
% hold on
% plot(mean(V(celltype.pool_stay & celltype.excit,:)))
% title('membrane voltage (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off
% subplot(2,1,2)
% plot(mean(V(celltype.pool_switch & celltype.inhib,:)))
% hold on
% plot(mean(V(celltype.pool_switch & celltype.excit,:)))
% title('membrane voltage (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% hold off

