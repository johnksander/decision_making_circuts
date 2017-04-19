clear
clc
close all
format compact

basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
addpath(fullfile(basedir,'helper_functions'))

%equal pools & currents (b set to 3) (WORKS!)

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
pool_options.p_conn = NaN; %randomized connection weights, not spatial connections
%--------------------------------------------------------------------------
%build the connectivity matrix
%make celltype logicals
celltype = celltype_logicals(pool_options);
%make these vectors into logical matricies
[EtoE,EtoI,ItoE] = connection_logicals(celltype,pool_options.num_cells);
%make connection scheme based off connection matricies
connection_scheme = EtoE | EtoI | ItoE;  %plan is:  EtoE | EtoI | ItoE;
%make synaptic connection matrix
W = rand(pool_options.num_cells); %use randomized connection weights
W(~connection_scheme) = 0; %filter out connections not allowed by scheme
%modify connection weights
W(EtoE) = .025 * W(EtoE);
W(ItoE) = 6 * W(ItoE);
W(EtoI) = .1 * W(EtoI);
%reorder weight matrix for column indexing in loop 
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
current_val = 80; %applied current !!NOTE!!: pulse must be added for different inhibitory & excitatory currents like he has it
noise_sigma = 400; %noise variance in picoamps
current_pulse = 'on'; %switch for adding transient pulses (also use for different excit & inhib currents)
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential, excitatory
Erev(celltype.inhib) = -62; %reversal potential, inhibitory
Gg = NaN(pool_options.num_cells,1); %max connductance
Gg(celltype.excit) = 1; %just setting these to 1, has no impact.
Gg(celltype.inhib) = 1; %just setting these to 1, has no impact.
Pr = NaN(pool_options.num_cells,1); %release probability
Pr(celltype.excit) = .01; %excitatory release probability
Pr(celltype.inhib) = .01; %inhibitory release probability
W(:,celltype.excit) = W(:,celltype.excit)/unique(Pr(celltype.excit)); %modify connection weights by Pr
Td = 5;%synaptic depression time constant, seconds
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50; %excitatory gating time constant, ms
Tsyn(celltype.inhib) = 5; %inhibitory gating time constant, ms
%----cell basics--------------
El = -62; %leak potential mV
Vreset = NaN(pool_options.num_cells,1); %reset potential mV
Vreset(celltype.excit) = -58; %excitatory
Vreset(celltype.inhib) = -56; %inhibitory
Gl = NaN(pool_options.num_cells,1); %leak conductance
Gl(celltype.excit) = 12; %excitatory
Gl(celltype.inhib) = 10; %inhibitory
Rm = 1./Gl; %resistance
Cm = 200; %cell capacity picofarads
spike_thresh = 20; %spike reset threshold (higher than Vth)
%---adaptation current--------
Vth = -50; %ALEIF spike threshold mV
delta_th = 2; %max voltage threshold, mV  (deltaVth in equation)
a = 2; %adaptation conductance, in nano Siemens
b = NaN(pool_options.num_cells,1); %increase adaptation current by max value, nano amps
b(celltype.excit & celltype.pool_stay) = 3; %excitatory- stay pool
%b(celltype.excit & celltype.pool_switch) = 0.5; %excitatory- switch pool
b(celltype.excit & celltype.pool_switch) = 3; %excitatory- switch pool
b(celltype.inhib) = 0; %inhibitory
Tsra = NaN(pool_options.num_cells,1);%adaptation current time constant, ms
Tsra(celltype.excit) = 300; %excitatory
Tsra(celltype.inhib) = 30; %inhibitory
%----timecourse---------------
tmax = 90000; %simulation end, what's the unit here?
timestep = .25; %what's the unit here?
timevec = 0:timestep:tmax;
switch current_pulse
    case 'on'
        diffEcurrent = 20; %he has different E & I currents, add this for that effect
        pulse_params = {celltype.excit,diffEcurrent,0,tmax};
end
%-------------------------------------------------------------------------
%preallocate variables

V = NaN(pool_options.num_cells,numel(timevec)); %membrane potential
V(:,1) = -80+40*rand(numel(V(:,1)),1); %inital value of membrane potential is random

Iapp = zeros(size(V)); %applied current vector
Iapp = Iapp + current_val;
noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment
Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise
switch current_pulse
    case 'on'
        Iapp = addpulse(pulse_params,Iapp,timevec); %add pulses
end

Isra = NaN(size(V)); %adaptation current
Isra(:,1) = 40*rand(numel(Isra(:,1)),1); %initialize at random

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
    D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
    Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
    
    if  sum(V(:,idx) > spike_thresh) > 0
        spiking_cells = V(:,idx) > spike_thresh;
        Isra(spiking_cells,idx) = Isra(spiking_cells,idx) + b(spiking_cells); %increase adaptation current by max value
        Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
            (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
        Sg(spiking_cells & celltype.inhib,idx) = 1; %synaptic gating for inhibitory cells 
        D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
        D(spiking_cells & celltype.inhib,idx) = 1; %synaptic depression for inhibitory cells 
        V(spiking_cells,idx) = Vreset(spiking_cells);
        spikes(spiking_cells,idx) = 1;
    end
    
    %[state,durations] = test4switch(Sg(:,idx),state,durations);
    
end

example_cells = [115 20];
timestep = timestep/1000; %resetting for figures below, if /1000 is even right?

figure(1)
spikeplot = make_spikeplot(spikes);
imagesc(spikeplot)
Xtick_seconds=0:.5:tmax/1000; %I think?
Xtick_pos = find(mod(timevec./1000,.5) == 0); %also-see above
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
title({'spikes','(spikes in matrix enlarged for visualization)'})
ylabel('cell')
xlabel('time (s)')
figure(2)
imagesc(Iapp)
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
title('current')
ylabel('cell')
xlabel('time (s)')

figure(3)
subplot(2,1,1)
plot(D(example_cells(1),:))
hold on
plot(D(example_cells(2),:))
hold off
title('synaptic depression (stay pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
subplot(2,1,2)
plot(D(example_cells(1) + 125,:))
hold on
plot(D(example_cells(2) + 125,:))
hold off
title('synaptic depression (switch pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))




figure(4)
subplot(2,1,1)
plot(Sg(example_cells(1),:))
hold on
plot(Sg(example_cells(2),:))
title('synaptic gating (stay pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
hold off
subplot(2,1,2)
plot(Sg(example_cells(1) + 125,:))
hold on
plot(Sg(example_cells(2) + 125,:))
title('synaptic gating (switch pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
hold off

figure(5)
subplot(2,1,1)
plot(V(example_cells(1),:))
hold on
plot(V(example_cells(2),:))
title('membrane voltage (stay pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
hold off
subplot(2,1,2)
plot(V(example_cells(1) + 125,:))
hold on
plot(V(example_cells(2) + 125,:))
title('membrane voltage (switch pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
hold off



figure(6)
subplot(2,1,1)
plot(Isra(example_cells(1),:))
hold on
plot(Isra(example_cells(2),:))
hold off
legend('inhibitory','excitatory')
title('Adaptation current (stay pool example cells)')
xlabel('time (s)')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
subplot(2,1,2)
plot(Isra(example_cells(1) + 125,:))
hold on
plot(Isra(example_cells(2) + 125,:))
hold off
legend('inhibitory','excitatory')
title('Adaptation current (switch pool example cells)')
xlabel('time (s)')
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))


figure(7)
subplot(2,1,1)
Irateplot = make_spikerate_plot(spikes,celltype.pool_stay & celltype.inhib,timestep);
Irateplot(isnan(Irateplot)) = 0;
plot(Irateplot,'linewidth',2)
hold on
Erateplot = make_spikerate_plot(spikes,celltype.pool_stay & celltype.excit,timestep);
Erateplot(isnan(Erateplot)) = 0;
plot(Erateplot,'linewidth',2)
title('stay pool cells')
legend('inhibitory','excitatory')
ylabel({'average cell rate (Hz)','(mean cell 1/ISI at each timepoint)'})
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
xlabel('time (s)')
hold off
subplot(2,1,2)
Irateplot = make_spikerate_plot(spikes,celltype.pool_switch & celltype.inhib,timestep);
Irateplot(isnan(Irateplot)) = 0;
plot(Irateplot,'linewidth',2)
hold on
Erateplot = make_spikerate_plot(spikes,celltype.pool_switch & celltype.excit,timestep);
Erateplot(isnan(Erateplot)) = 0;
plot(Erateplot,'linewidth',2)
title('switch pool cells')
legend('inhibitory','excitatory')
ylabel({'average cell rate (Hz)','(mean cell 1/ISI at each timepoint)'})
set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
xlabel('time (s)')
hold off


figure(8)
histogram(Iapp(1,:))
title('current distribution, cell #1')




% disp('---during pulse:')
% pulse_time = timevec >= ptime(1) & timevec <= ptime(2);
% Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
% Erate = mean(Erate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
% Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
% Irate = mean(Irate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
disp('---during t+ two seconds (stay):')
pulse_time = timevec >= 2;
Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
Erate = mean(Erate) / (sum(pulse_time) * timestep);
disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
Irate = mean(Irate) / (sum(pulse_time) * timestep);
disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
disp('---during t+ two seconds (switch):')
pulse_time = timevec >= 2;
Erate = sum(spikes(celltype.pool_switch & celltype.excit,pulse_time),2);
Erate = mean(Erate) / (sum(pulse_time) * timestep);
disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
Irate = sum(spikes(celltype.pool_switch & celltype.inhib,pulse_time),2);
Irate = mean(Irate) / (sum(pulse_time) * timestep);
disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))



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

