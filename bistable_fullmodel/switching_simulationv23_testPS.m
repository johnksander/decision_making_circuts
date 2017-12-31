clear
clc
close all
format compact
rng('shuffle') %this is probably important...

basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project';
addpath(fullfile(basedir,'helper_functions'))
resfile = 'PS_parsweep_baseline_5';
resdir = 'parsweep_baseline';
options = load(fullfile(basedir,'Results',resdir,resfile,resfile));
options = options.options;

%parameterization from v22.
%the network should satisfy:
%   a) slow/no switching when all things are equal
%   b) stimulus casues fast switching
%very well.

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 250;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.8 .2]; %proportion excitable % inhibitory
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
W(W > 0 & EtoE) = options.EtoE;
W(W > 0 & ItoE) = options.ItoE;
W(W > 0 & EtoI) = options.EtoI;
%reorder weight matrix for column indexing in loop
W = reorder_weightmat(W,celltype);
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential, excitatory
Erev(celltype.inhib) = -70e-3; %reversal potential, inhibitory
Gg = NaN(pool_options.num_cells,1); %max connductance
Gg(celltype.excit) = 10e-9; %excitatory max connductance microSiemens
Gg(celltype.inhib) = 10e-9; %inhibitory max connductance microSiemens
Pr = NaN(pool_options.num_cells,1); %release probability
Pr(celltype.excit) = .2; %excitatory release probability
Pr(celltype.inhib) = .2; %inhibitory release probability
Td = .3;%synaptic depression time constant, seconds
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50e-3; %excitatory gating time constant, ms
Tsyn(celltype.inhib) = 10e-3; %inhibitory gating time constant, ms
%----cell basics--------------
El = -70e-3; %leak potential mV
Ek = -80e-3; %potassium potential mV
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%----noisy input--------------
Tau_ext = NaN(pool_options.num_cells,1); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 2e-3;
Tau_ext(celltype.inhib) = 5e-3;
initGext = 10e-9; %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9; %increase noisy conducrance, nano Siemens
Rext = 1400; %poisson spike train rate for noise, Hz
%Rstim = 30; %rate for stimulus input spikes
Rstim = 0; %rate for stimulus input spikes
stim_targs = celltype.excit & celltype.pool_switch;
%---adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation current time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----timecourse---------------
tmax = 100; %simulation end (s)
timestep = .25e-3; %.01 milisecond timestep
timevec = 0:timestep:tmax;
%-------------------------------------------------------------------------
%preallocate variables

V = NaN(pool_options.num_cells,numel(timevec)); %membrane potential
V(:,1) = El; %inital value of membrane potential is leak potential

Lext = Rext * timestep; %poisson lambda for noisy conductance

Gext = NaN(size(V)); %noisy conductance
Gext(:,1) = initGext; %initialize at leak conductance

Gsra = NaN(size(V)); %adaptation conductance
Gsra(:,1) = 0; %initialize at zero

Sg = NaN(size(V)); %synaptic gating
Sg(:,1) = 0; %initalize at zero??
D = NaN(size(V)); %synaptic depression
D(:,1) = 1; %initalize at one

spikes = zeros(size(V));

%---state tracker-------------
state = struct();
state.stay = logical([1 0]);
state.switch = logical([0 1]);
durations = cell(1,2); %record duration times (stay, switch)
state.now = state.stay; %pick one state to start with, add pulse to that pool to be sure (make sure this is consistent)
state.stim_labels = {'A','B'};
state.current_stimulus = logical([1 0]); %initialize in stim A
state.count = 0;
state.pools2compare = [celltype.pool_stay & celltype.excit,...
    celltype.pool_switch & celltype.excit]; %pass in this format, avoid many computations
state.ready_mintime = .4 / timestep; %minimum time for ready2go check
if options.force_back2stay
    state.back2stay_min = .5 / timestep; %wait in a stay state for this long before forcing switch back
end
%---stimuli info--------------
stim_info = struct();
stim_info.targ_cells = stim_targs; %Estay
stim_info.stimA_lambda = Rstim * timestep; %poisson lambda for stimulus conductance
stim_info.stimB_lambda = Rstim * timestep;
stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well

timepoint_idx = 1;

for idx = 2:numel(timevec)
    
    timepoint_idx = timepoint_idx + 1;
    
    I = (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
    I = I + (unique(Erev(celltype.inhib)) - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*unique(Gg(celltype.inhib));
    dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm)...
        + (Gsra(:,idx-1).*(Ek-V(:,idx-1)))...
        + (Gext(:,idx-1).*(Erev-V(:,idx-1))) + I;
    
    V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
    
    Gsra(:,idx) = Gsra(:,idx-1) - ((Gsra(:,idx-1)./Tsra) .* timestep); %adaptation conductance
    Gext(:,idx) = Gext(:,idx-1) - ((Gext(:,idx-1)./Tau_ext) .* timestep); %noisy conductance
    D(:,idx) = D(:,idx-1) + (((1 - D(:,idx-1))./Td) .* timestep); %synaptic depression
    Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
    
    if  sum(V(:,idx) > spike_thresh) > 0
        spiking_cells = V(:,idx) > spike_thresh;
        Gsra(spiking_cells,idx) = Gsra(spiking_cells,idx) + detlaGsra; %adaptation conductance
        Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
            (Pr(spiking_cells).*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
        D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr(spiking_cells)); %synaptic depression
        V(spiking_cells,idx) = Vreset;
        spikes(spiking_cells,idx) = 1;
    end
    
    [state,durations] = test4switch(Sg(:,idx),state,durations);
    
    
    ext_spikes = poissrnd(Lext,pool_options.num_cells,1); %external spikes to noise conductance
    stim_spikes = timepoint_stimulus(stim_info,state);
    ext_spikes = ext_spikes + stim_spikes;
    if options.force_back2stay
        ext_spikes = back2stay(ext_spikes,state,celltype); %if in switch state, force back to a stay state
    end

    
    Gext(:,idx) = Gext(:,idx) + (deltaGext.*ext_spikes); %don't have to index, they get an increase or zero
    
end

example_cells = [60 20];

figure(1)
spikeplot = make_spikeplot(spikes);
imagesc(spikeplot)

Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title({'spikes','(spikes in matrix enlarged for visualization)'})
ylabel('cell')
xlabel('time (s)')
% figure(2)
% imagesc(Iapp)
% title('current')
% ylabel('cell')
% xlabel('time (s)')
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);


% figure(3)
% subplot(2,1,1)
% plot(D(example_cells(1),:))
% hold on
% plot(D(example_cells(2),:))
% hold off
% title('synaptic depression (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% subplot(2,1,2)
% plot(D(example_cells(1) + sum(celltype.pool_stay),:))
% hold on
% plot(D(example_cells(2) + sum(celltype.pool_stay),:))
% hold off
% title('synaptic depression (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);



% figure(4)
% subplot(2,1,1)
% plot(Sg(example_cells(1),:))
% hold on
% plot(Sg(example_cells(2),:))
% title('synaptic gating (stay pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% hold off
% subplot(2,1,2)
% plot(Sg(example_cells(1) + sum(celltype.pool_stay),:))
% hold on
% plot(Sg(example_cells(2) + sum(celltype.pool_stay),:))
% title('synaptic gating (switch pool example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% hold off

figure(5)
subplot(2,1,1)
plot(V(example_cells(1),:))
hold on
plot(V(example_cells(2),:))
title('membrane voltage (stay pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
hold off
subplot(2,1,2)
plot(V(example_cells(1) + sum(celltype.pool_stay),:))
hold on
plot(V(example_cells(2) + sum(celltype.pool_stay),:))
title('membrane voltage (switch pool example cells)')
xlabel('time (s)')
legend('inhibitory','excitatory')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
hold off



% figure(6)
% subplot(2,1,1)
% plot(Gsra(example_cells(1),:))
% hold on
% plot(Gsra(example_cells(2),:))
% hold off
% legend('inhibitory','excitatory')
% title('Adaptation conductance (stay pool example cells)')
% xlabel('time (s)')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% subplot(2,1,2)
% plot(Gsra(example_cells(1) + sum(celltype.pool_stay),:))
% hold on
% plot(Gsra(example_cells(2) + sum(celltype.pool_stay),:))
% hold off
% legend('inhibitory','excitatory')
% title('Adaptation conductance (switch pool example cells)')
% xlabel('time (s)')
% Xticks = num2cell(get(gca,'Xtick'));
% Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);

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
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
xlabel('time (s)')
set(gca,'Xlim',[0 numel(Erateplot)])
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
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
xlabel('time (s)')
set(gca,'Xlim',[0 numel(Erateplot)])
hold off

% figure(8)
% histogram(Iapp(1,:))
% title('current distribution, cell #1')

figure(9)
subplot(2,1,1)
plot(Gext(example_cells(1),:))
hold on
plot(Gext(example_cells(2),:))
hold off
legend('inhibitory','excitatory')
title('Noisy background conductance (stay pool example cells)')
xlabel('time (s)')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
subplot(2,1,2)
plot(Gext(example_cells(1) + sum(celltype.pool_stay),:))
hold on
plot(Gext(example_cells(2) + sum(celltype.pool_stay),:))
hold off
legend('inhibitory','excitatory')
title('Noisy background conductance (switch pool example cells)')
xlabel('time (s)')
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%i',round(x*timestep)),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);




% disp('---during pulse:')
% pulse_time = timevec >= ptime(1) & timevec <= ptime(2);
% Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
% Erate = mean(Erate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
% Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
% Irate = mean(Irate) / (sum(pulse_time) * timestep);
% disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
cuttoff_time = 18;
disp(['---during t+ ' num2str(cuttoff_time) ' seconds (stay):'])
pulse_time = timevec >= cuttoff_time;
Erate = sum(spikes(celltype.pool_stay & celltype.excit,pulse_time),2);
Erate = mean(Erate) / (sum(pulse_time) * timestep);
disp(sprintf('mean excitatory spike rate %.1f Hz',Erate))
Irate = sum(spikes(celltype.pool_stay & celltype.inhib,pulse_time),2);
Irate = mean(Irate) / (sum(pulse_time) * timestep);
disp(sprintf('mean inhibitory spike rate %.1f Hz',Irate))
disp(['---during t+ ' num2str(cuttoff_time) ' seconds (switch):'])
pulse_time = timevec >= cuttoff_time;
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

