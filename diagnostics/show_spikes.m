

pool_options.num_cells = 1;
tmax = .1;
Rext = 15e3;
timestep = .02e-3;
num_timepoints = round(tmax / timestep) + 1; %same as numel(0:dt:Tmax)
Lext = Rext * timestep; %poisson lambda for noisy conductance

celltype.excit = true;
celltype.inhib = false;


W = 0;
W_Ex = W(:,celltype.excit);
W_In = W(:,celltype.inhib);

%set up simulation parameters
%--------------------------------------------------------------------------
%----cell basics---------------
Erev = 0; %reversal potential, excitatory
Irev = -70e-3; %reversal potential, inhibitory
Gg = 10e-9; %max conductance microSiemens
Tsyn_Ex = 50e-3; %excitatory gating time constant, ms
Tsyn_In = 10e-3; %inhibitory gating time constant, ms
El = -70e-3; %leak potential mV
Ek = -80e-3; %potassium potential mV
Vreset = -80e-3; %reset potential mV
Rm = 100e6; %resistance megaohms
Gl = 1/Rm; %leak conductance
Cm = 100e-12; %cell capacity picofarads
spike_thresh = 20e-3; %spike reset threshold (higher than Vth)
%----noisy input---------------
Tau_ext = NaN(pool_options.num_cells,1); %noisy conductance time constant, ms
Tau_ext(celltype.excit) = 3.5e-3; % was 2e-3
Tau_ext(celltype.inhib) = 2e-3; % was 5e-3
initGext = 10e-9; %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9; %increase noisy conducrance, nano Siemens
Rext = 1540; %poisson spike train rate for noise, Hz
%----adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----depression----------------
Pr = .1; %release probability
Td.fast = .3;%synaptic depression time constant, seconds (fast)
Td.slow = 10;%slow time constant
fD = .05; %availability factor: (max docked vessicles / max pooled vessicles)

%---membrane potential--------
V = zeros(pool_options.num_cells,1);
V = V + El;%inital value of membrane potential is leak potential
%---stimuli info--------------
%for trials, work out passing trial_stimuli indexed by trialidx to init_statevar()
%---noisy conductance---------
ext_inds.I = 1;
ext_inds.E = 2;
Gext_Ex = initGext + zeros(pool_options.num_cells,1); %initialize at leak conductance
Gext_In = initGext + zeros(pool_options.num_cells,1);
%---adaptation conductance----
Gsra = zeros(size(V));
%---gating & depression-------
Sg_Ex = zeros(sum(celltype.excit),1);
Sg_In = zeros(sum(celltype.inhib),1);
Dfast = ones(size(V)); %synaptic depression: fast
Dslow = ones(size(V)); %synaptic depression: slow
%---spikes--------------------
timepoint_counter = 1;

Vrec = NaN(pool_options.num_cells,num_timepoints-1);
Nspike = 0;

while timepoint_counter < num_timepoints
    
    timepoint_counter = timepoint_counter+1;
    
    WS_Ex = W_Ex*Sg_Ex;
    WS_In = W_In*Sg_In;
    
    E_diff = Erev-V;
    I_diff = Irev-V;
    I = E_diff.*WS_Ex.*Gg + I_diff.*WS_In.*Gg + Gext_Ex.*E_diff + Gext_In.*I_diff;
    dVdt = ((El - V + (delta_th.*exp((V - Vth)./delta_th)))./Rm) + (Gsra.*(Ek - V)) + I;
    V = ((dVdt./Cm) .* timestep) + V;
    
    Gext_Ex = Gext_Ex - (Gext_Ex./Tau_ext) .* timestep;%noisy conductance
    Gext_In = Gext_In - (Gext_In./Tau_ext) .* timestep;
    Gsra = Gsra - (Gsra./Tsra) .* timestep;%adaptation conductance
    Sg_Ex = Sg_Ex - (Sg_Ex./Tsyn_Ex) .* timestep;%synaptic gating
    Sg_In = Sg_In - (Sg_In./Tsyn_In) .* timestep;
    Ddiff = Dslow - Dfast;
    Dfast = Dfast + ((Ddiff ./ Td.fast) .* timestep);%fast syn. depression
    Dslow = Dslow + timestep .* ( ((1 - Dslow) ./ Td.slow) ...
        - fD.* (Ddiff ./ Td.fast)  ); %slow vessicle replacement
    
    spiking_cells = V > spike_thresh;
    if sum(spiking_cells) > 0
        Gsra(spiking_cells) = Gsra(spiking_cells) + detlaGsra; %adaptation conductance
        %synaptic gating, updates with Dfast vessicles
        Espike = spiking_cells(celltype.excit);
        Ispike = spiking_cells(celltype.inhib);
        Sg_Ex(Espike) = Sg_Ex(Espike) + (Pr.* Dfast(spiking_cells & celltype.excit).*(1-Sg_Ex(Espike)));
        Sg_In(Ispike) = Sg_In(Ispike) + (Pr.* Dfast(spiking_cells & celltype.inhib).*(1-Sg_In(Ispike)));
        %depression update (docked vessicles released)
        Dfast(spiking_cells) = Dfast(spiking_cells) - (Pr.* Dfast(spiking_cells));
        V(spiking_cells) = Vreset;
        Nspike = Nspike + 1;
    end
    
    
    
    
    %input spikes: noise & stimulus
    %---noisy spiking input from elsewhere
    ext_spikes = poissrnd(Lext,pool_options.num_cells,2);
    
    
    
    %update Gexternal. Don't have to index, they get an increase or zero
    Gext_Ex = Gext_Ex + deltaGext.*ext_spikes(:,ext_inds.E);
    Gext_In = Gext_In + deltaGext.*ext_spikes(:,ext_inds.I);
    
    
    Vrec(:,timepoint_counter) = V;
    
end

lnsz = 10;
lncol = 'k';

Hz = Nspike / tmax;
plot(Vrec,'LineWidth',lnsz,'Color',lncol); 
pl_info = sprintf('%i Hz for %.1fs',round(Hz),tmax);
fprintf('%s\n',pl_info)
fn = strrep(pl_info,'.','');
fn = strrep(fn,' ','_');
axis tight
set(gca,'TickLength',[0,0])
set(gca,'XTickLabel',{})
set(gca,'YTickLabel',{})
set(gcf,'Renderer','painters')
print(fn,'-dpng','-r600');
ax = gca;
exportgraphics(gca,[fn,'.png'],'Resolution',300) 
