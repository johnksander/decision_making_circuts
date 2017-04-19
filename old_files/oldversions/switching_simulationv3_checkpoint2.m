clear
clc
format compact

%This works only with E-E connections (within pool excitatory feedback)
%and no cross pool connections. Also, only with Tsyn (synaptic gating time
%constant) set to 10e-3 


basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
addpath(fullfile(basedir,'helper_functions'))

num_cells = 160;
pool_options.sz_pools = [.5 .5]; %proportion stay & switch
pool_options.sz_EI = [.5 .5]; %proportion excitable % inhibitory 
[pool_stay,pool_switch,excit,inhib] = celltype_logicals(num_cells,pool_options); %make celltype logicals 


p_conn = .5; %connection probability 50%
W = rand(num_cells) < p_conn; %connection matrix
%I remember something about adding noise here?

%make these vectors into logical matricies
%pool_stay_mat = repmat(pool_stay,1,num_cells) | repmat(pool_stay',num_cells,1); %maybe figure out how these need to be used later
%pool_switch_mat = repmat(pool_switch,1,num_cells) | repmat(pool_switch',num_cells,1);
EtoE =  repmat(excit,1,num_cells) & repmat(excit',num_cells,1);
EtoI =  repmat(excit,1,num_cells) & repmat(inhib',num_cells,1);
ItoE =  repmat(inhib,1,num_cells) & repmat(excit',num_cells,1);
self_connection = repmat(pool_stay,1,num_cells) == repmat(pool_stay',num_cells,1);
EtoE = EtoE & self_connection; %only self-excitation within group
EtoI = EtoI & ~self_connection; %E to I connection is across group
ItoE = ItoE & self_connection; %I to E connection is within group

%connection_scheme = EtoE | EtoI | ItoE;
connection_scheme = EtoE ;
%filter out connections not allowed by scheme
W = double(W & connection_scheme);
%increase depression
%W(W > 0 & ItoE) = 10; %watch out for how this is selected (W ~=0) if you add noise or something
%W(W > 0 & EtoE) = 1; %watch out for how this is selected (W ~=0) if you add noise or something


%set up parameters
%--------------------------------------------------------------------------
%----cell connections---------
Erev = NaN(num_cells,1); %reversal potential vector
Erev(excit) = -70e-3; %reversal potential from inhibitory synapses (i.e. inhibitory synapse rev is -70mV)
Erev(inhib) = 0; %reversal potential from excitatory synapses
Gg = 1e-6; %max connductance microSiemens
Pr = .2; %release probability
Td = .2;%synaptic depression time constant, seconds
Tsyn = NaN(num_cells,1); %gating time constant vector 
Tsyn(excit) = 50e-3; %excitatory gating time constant, ms (applies to all excitatory synapses!)
Tsyn(excit) = 10e-3; %excitatory gating time constant, ms 
Tsyn(inhib) = 10e-3; %excitatory gating time constant, ms 
%----cell basics--------------
El = -70e-3; %leak potential mV
Vth = -50e-3; %spike threshold mV
Vreset = -80e-3; %reset potential mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
%---adaptation current--------
a = 3e-9; %adaptation conductance, in nano Siemens
b = .02e-9;%increase adaptation current by max value, nano amps
Tsra = 200e-3; %adaptation current time constant, ms
%----timecourse---------------
tmax = 3; %simulation end (s)
timestep = .25e-3; %.01 milisecond timestep
timevec = 0:timestep:tmax;
%----current and noise--------
current_val = 1e-9; %applied current in nano amps
%current_val =  0
noise_sigma = 50e-12; %noise variance in picoamps
%-------------------------------------------------------------------------
%preallocate variables
V = NaN(num_cells,numel(timevec)); %membrane potential
V(:,1) = El; %inital value of membrane potential is leak potential

Iapp = zeros(size(V)); %applied current vector
Iapp = Iapp + current_val;
noise_sigma = noise_sigma/sqrt(timestep); %take care of timestep adjustment
Iapp = Iapp + (randn(size(Iapp)).*noise_sigma); %add noise

Isra = NaN(size(V)); %adaptation current
Isra(:,1) = 0; %initialize at zero??

Sg = NaN(size(V)); %synaptic gating
Sg(:,1) = 0; %initalize at zero??
D = NaN(size(V)); %synaptic depression
D(:,1) = 1; %initalize at one

spikes = zeros(size(V));


for idx = 2:numel(timevec)
    
    
    if timevec(idx) > .75 & timevec(idx) < 1.25 %current pulse to switch pool
        Iapp(pool_switch,idx) = Iapp(pool_switch,idx) + 5e-9; 
    end
    
    Isyn = W'*Sg(:,idx-1).*Gg;
    %all together now
    dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm)...
        + (Isyn.*(Erev - V(:,idx-1))) - Isra(:,idx-1) + Iapp(:,idx-1); %noise is already in Iapp
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
        D(spiking_cells,idx) = D(spiking_cells,idx).*(1-Pr); %synaptic depression
        V(spiking_cells,idx) = Vreset;
        spikes(spiking_cells,idx) = 1;
    end
end

imagesc(spikes)

