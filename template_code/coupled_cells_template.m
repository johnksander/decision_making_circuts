clear
clc
format compact

%1)
%set up params
Cm = 1e-9; %cell capacity nanofarads
Rm = 10e6; %resistance megaohms
E = -70e-3; %leak potential miliVolts (find out how this is different from Erev)
Vth = -54e-3; %spike threshold mV
Vr = -80e-3; %reset potential mV
Erev = -70e-3; %reverse potential miliVolts
G = 1e-6; %max connductance microSiemens
Tsyn = 10e-3; %gating time constant, ms 
Td = .2;%synaptic depression time constant, seconds
baselineIapp = 2e-9; %baseline applied current, nano amps
Pr = .2; 
%Pr = .25; 
%I can't get these to switch with a current pulse of 1nA unless Pr = .25(instead of 1 in the question)


%2)
timestep = .01e-3; %.01 milisecond timestep
timevec = 0:timestep:6; %tmax is 6 seconds
Iapp = zeros(2,numel(timevec)); %applied current vector, one row per cell
Iapp = Iapp + baselineIapp; %apply baseline current to both
%current_pulse = 1e-9; %additional current pulse, nano amps
current_pulse = 3e-9; %additional current pulse, nano amps
%I can't get these to switch unless I up the current pulse to 3nA (instead of 1 in the question)
pulse_time = timevec <= 100e-3; %current pulse for first 100 ms
Iapp(1,pulse_time) = Iapp(1,pulse_time) + current_pulse; %apply to first cell
pulse_time = timevec >= 3 & timevec <= 3.1; %current pulse for 100 ms starting at 3 seconds
Iapp(2,pulse_time) = Iapp(2,pulse_time) + current_pulse; %apply to second cell

V = NaN(size(Iapp)); %membrane potential, both cells
V(:,1) = E; %inital value of membrane potential is leak potential
spikes = zeros(size(V));
%S is synaptic gating
%D is synaptic depression'
D = NaN(size(V));
Sg = NaN(size(V));
Sg(:,1) = 0; %initalize at zero??
noise_sigma = 0; %no noise for this run
noise = randn(size(timevec)) * noise_sigma ./ sqrt(timestep); %noise is the same for both cells? (also, the * or / thing from previous chapt)
swap_sg = [2;1]; %for putting s1 in the cell 2 membrane potential equation, etc

for idx = 2:numel(timevec)
    
    D(:,idx) = D(:,idx-1) - ((D(:,idx-1)./Td) .* timestep); %synaptic depression
    Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
    dVdt = ((E - V(:,idx-1))./Rm) + (G.*Sg(swap_sg,idx-1).*(Erev - V(:,idx-1))) + Iapp(:,idx-1) + noise(idx-1);
    dVdt = dVdt ./ Cm; %membrane potential rate of change
    V(:,idx) = (dVdt .* timestep) + V(:,idx - 1);
    %G, Erev, and D are set as a constants for both cells here
    %noise is also the same for both cells?
    
    if sum(V(:,idx) > Vth) > 0
        spiking_cells = V(:,idx) > Vth;
        Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
            (Pr.*D(spiking_cells,idx).*(1-Sg(spiking_cells,idx))); %synaptic gating
        D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr); %synaptic depression 
        %dont reset D for this run
        V(spiking_cells,idx) = Vr;
        spikes(spiking_cells,idx) = 1;
    end
end

figure(1)
plot(V(1,:),'linewidth',2)
hold on
plot(V(2,:),'linewidth',2)
hold off
legend('Cell #1','Cell #2')
ylabel('Vm')
xlabel('time')

figure(2)
plot(Sg(1,:),'linewidth',2)
hold on
plot(Sg(2,:),'linewidth',2)
hold off
legend('Cell #1','Cell #2')
ylabel('Synaptic gating')
xlabel('time')

figure(3)
plot(spikes(1,:),'linewidth',2)
hold on
plot(spikes(2,:),'linewidth',2)
hold off
legend('Cell #1','Cell #2')
ylabel('spikes y/n')
xlabel('time')











