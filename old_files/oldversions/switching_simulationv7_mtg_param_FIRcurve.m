clear
clc
close all
format compact

basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
addpath(fullfile(basedir,'helper_functions'))

%set up the circut
%--------------------------------------------------------------------------
%----circit parameters--------
pool_options.num_cells = 200;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.5 .5]; %proportion excitable % inhibitory
pool_options.p_conn = .25; %connection probability 50%
%--------------------------------------------------------------------------
%build the connectivity matrix
%make celltype logicals
celltype = celltype_logicals(pool_options);
%make these vectors into logical matricies
[EtoE,EtoI,ItoE] = connection_logicals(celltype,pool_options.num_cells);
%make connection scheme based off connection matricies
%connection_scheme = EtoE | EtoI | ItoE;  %plan is:  EtoE | EtoI | ItoE;
connection_scheme = EtoE;
%make synaptic connection matrix
W = rand(pool_options.num_cells) < pool_options.p_conn; %connection matrix
W = double(W & connection_scheme); %filter out connections not allowed by scheme
%modify connection weights
W(W > 0 & EtoE) = .0;
% W(W > 0 & ItoE) = .2; 
% W(W > 0 & EtoI) = .05; 
%--------------------------------------------------------------------------

%set up simulation parameters
%--------------------------------------------------------------------------
%----current and noise--------
%current_val = 1e-9; %applied current in nano amps
%noise_sigma = 50e-12; %noise variance in picoamps
current_val = 0; %applied current in nano amps
noise_sigma = 0; %noise variance in picoamps
%----cell connections---------
Erev = NaN(pool_options.num_cells,1); %reversal potential vector
Erev(celltype.excit) = 0; %reversal potential from inhibitory synapses (i.e. inhibitory synapse rev is -70mV)
Erev(celltype.inhib) = -70e-3; %reversal potential from excitatory synapses
Gg = NaN(pool_options.num_cells,1); %gating time constant vector
Gg(celltype.excit) = .5e-6; %excitatory max connductance microSiemens
Gg(celltype.inhib) = .5e-6; %inhibitory max connductance microSiemens
Pr = .2; %release probability
Td = .2;%synaptic depression time constant, seconds
Tsyn = NaN(pool_options.num_cells,1); %gating time constant vector
Tsyn(celltype.excit) = 50e-3; %excitatory gating time constant, ms 
Tsyn(celltype.inhib) = 10e-3; %inhibitory gating time constant, ms
%----cell basics--------------
El = -70e-3; %leak potential mV
Vth = -50e-3; %spike threshold mV
Vreset = -80e-3; %reset potential mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
%---adaptation current--------
a = 2e-9; %adaptation conductance, in nano Siemens
b = .02e-9;%increase adaptation current by max value, nano amps 
Tsra = 200e-3; %adaptation current time constant, ms
%----timecourse---------------
tmax = 5; %simulation end (s)
timestep = .25e-3; %.01 milisecond timestep
timevec = 0:timestep:tmax;
% pulse_current = 20e-9;
% pulse_params = {celltype.pool_stay,pulse_current,1,1.25;...
%     celltype.pool_stay,pulse_current,3,5}; %logical, current, start, end
% pulse_params = {celltype.pool_stay,pulse_current,1,4};...

%-------------------------------------------------------------------------


FIRcurrents = .1e-9:.05e-9:.6e-9;
num_trials = numel(FIRcurrents);
spikerates = NaN(pool_options.num_cells,num_trials); %recorded in Hz

for trialidx = 1:num_trials
    %preallocate variables
    
    pulse_current = FIRcurrents(trialidx);
    pulse_params = {celltype.pool_stay,pulse_current,0,tmax}; %give pulse to whole timecourse
    
    V = NaN(pool_options.num_cells,numel(timevec)); %membrane potential
    V(:,1) = El; %inital value of membrane potential is leak potential
    
    Iapp = zeros(size(V)); %applied current vector
    Iapp = Iapp + current_val;
    noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment
    Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise
    Iapp = addpulse(pulse_params,Iapp,timevec); %add pulses
    
    Isra = NaN(size(V)); %adaptation current
    Isra(:,1) = 0; %initialize at zero??
    
    Sg = NaN(size(V)); %synaptic gating
    Sg(:,1) = 0; %initalize at zero??
    D = NaN(size(V)); %synaptic depression
    D(:,1) = 1; %initalize at one
    
    spikes = zeros(size(V));
    
    
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
        
        if  sum(V(:,idx) > Vth) > 0
            spiking_cells = V(:,idx) > Vth;
            Isra(spiking_cells,idx) = Isra(spiking_cells,idx) + b; %increase adaptation current by max value
            Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
                (Pr.*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
            %can add inhbibitory sg reset to 1
            %can stop super high firing states, quick spikes don't do anything
            D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr); %synaptic depression
            V(spiking_cells,idx) = Vreset;
            spikes(spiking_cells,idx) = 1;
        end
        
    end
    spikerates(:,trialidx) = sum(spikes,2) ./ tmax; %recorded in Hz
    disp(sprintf('---Trial #%i complete',trialidx))
end


poolrates = spikerates(celltype.pool_stay,:);
Erates = mean(spikerates(celltype.pool_stay & celltype.excit,:));
Irates = mean(spikerates(celltype.pool_stay & celltype.inhib,:));
plot(FIRcurrents,Erates,'linewidth',2)
hold on
plot(FIRcurrents,Irates,'linewidth',2)
hold off
title({'FIR curve','Only one pool with no connections- mtg params'})
ylabel('average cell spike Rate (Hz)')
Xcurrents = FIRcurrents / 1e-9;
Xcurrents = Xcurrents(1:5:end);
set(gca,'XTickLabel', arrayfun(@num2str, Xcurrents(:), 'UniformOutput', false))
xlabel('Iapp (nA)')
legend('Excitatory cells','Inhibitory cells','location','southeast')
print(fullfile(basedir,'switching_simv7_mtg_params_FIRcurve'),'-djpeg')
%I guess try taking a single neuron or something...



%
% figure(1)
% imagesc(spikes)
% Xtick_seconds=0:.5:tmax;
% Xtick_pos = find(mod(timevec,.5) == 0);
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% title('spikes')
% ylabel('cell')
% xlabel('time (s)')
% figure(2)
% imagesc(Iapp)
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% title('current')
% ylabel('cell')
% xlabel('time (s)')
% figure(3)
% plot(D(60,:))
% hold on
% plot(D(20,:))
% hold off
% title('synaptic depression (example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))
% figure(4)
% plot(Sg(60,:))
% hold on
% plot(Sg(20,:))
% title('synaptic gating (example cells)')
% xlabel('time (s)')
% legend('inhibitory','excitatory')
% set(gca,'Xdir','normal','XTick', Xtick_pos,'XTickLabel', arrayfun(@num2str, Xtick_seconds(:), 'UniformOutput', false))

% figure(2)
% histogram(durations{1} * timestep)
% hold on
% histogram(durations{2} * timestep)
% legend('stay','switch')
% title('state durations')
% ylabel('frequency')
% xlabel('time (s)')
% hold off



