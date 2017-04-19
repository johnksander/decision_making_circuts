clear
clc
format compact


%set up params
El = -75e-3; %leak potential mV
Vth = -50e-3; %spike threshold mV
Vreset = -80e-3; %reset potential mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
a = 2e-9; %adaptation conductance, in nano Siemens
b = .02e-9;%adaptation current, nano amps
Tsra = 200e-3; %adaptation current time constant, ms
tmax = 3; %simulation end (s) 
timestep = .25e-3; %.01 milisecond timestep
timevec = 0:timestep:tmax; %tmax is 1.5 seconds
current_val = 3e-9; %applied current in nano amps
Iapp = zeros(size(timevec)); %applied current vector
Iapp = Iapp + current_val;
noise_sigma = 400/sqrt(timestep); %current noise (incorproates timestep here)

keyboard



V = NaN(size(timevec)); %membrane potential
V(1) = El; %inital value of membrane potential is leak potential
Isra = NaN(size(timevec)); %adaptation current
Isra(1) = 0; %initialize at zero??
spikes = zeros(size(V));
for idx = 2:numel(timevec)
    dt_Isra = ((a * (V(idx-1) - El)) - Isra(idx-1)) / Tsra ; %adaptation current rate of change
    Isra(idx) = (dt_Isra * timestep) + Isra(idx-1);
    dVdt = (Gl * (El - V(idx-1) + (delta_th * exp((V(idx-1) - Vth)/delta_th)))) - Isra(idx-1) + Iapp(idx-1);
    dVdt = dVdt / Cm; %membrane potential rate of change
    V(idx) = (dVdt * timestep) + V(idx - 1);
    if V(idx) > Vth
        V(idx) = Vreset;
        Isra(idx) = Isra(idx) + b; %increase adaptation current by max value
        spikes(idx) = 1;
    end
end

subplot(2,1,1)
plot(Iapp,'linewidth',2)
ylabel('current')
xlabel('time')
subplot(2,1,2)
plot(V,'linewidth',2)
ylabel('Vm')
xlabel('time')
print('q2a_figure','-djpeg')






