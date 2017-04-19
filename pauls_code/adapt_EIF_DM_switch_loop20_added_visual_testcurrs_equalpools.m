% adapting_DM.m
% decision-making network with adaptation to implement switching
clear
clc
format compact



basefile = 'DM_switch_STP_20';

duration_file = strcat(basefile, 'dur');


NE=200; %excitatory neurons
NI=50; %inhibitory neurons
WEE0 = 5/NE;
WEI0 = 20/NE; %stronger e- i than e-e
WIE0 = 300/NI;%inhibitory - excitatory strongest %*6
WII0 = 0/NI; %inhibitory - inhibitory - no connections
%he has a mismatched number of ecitatiory & inhibitory neurons

% WEE0 = 0/NE;
% WEI0 = 0/NE;
% WIE0 = 0/NI;
% WII0 = 0/NI;

aE= 2; %time scale of recovery (typical = .02, so this is really fast)
aI = 2; %time scale of recovery (very fast) <- these are for adaptation current
CE = 200; %cell capacity excitatory
CI = 200; %cell capacity inhibitory
gLE = 12; %leak conductance excitatory
gLI = 10; %leak conductance inhibitory
ELE = -62; %leak potential excitatory
ELI = -62; %leak potential inhibitory
VTE = -50; % ALEIF Vth
VTI = -50;
deltaTE = 2; %ALEIF delta-th parameter
deltaTI = 2;
tauwE = 300; %tauSRA adaptation conductance time constant
tauwI = 30;

%Seems like these are the only things unbalancing the pools 
%bE1 = 3; %b looks like it might be Isra + b after spike
%bE2 = 0.5;
bE1 = 3;
bE2 = bE1;

bI = 0;
VRE = -58;
VRI = -56;



%he had it at 900000, but my computer has real life memory...
tmax = 90000;

dt = .25;   % in ms
tvector = dt:dt:tmax;
disable_init_pulse = 'yes';

switch disable_init_pulse
    case 'yes'
        t_init = 0; %disable pulse
        disp('initial pulse disabled')
    case 'no'
        t_init = 200; %some kind of pulse timing
end

Nwindow = floor(10000/dt); %unsure

I0 = 100; %unsure
I0I = 80;
Iswitch0 = 120;

Isigma = 400/sqrt(dt); %current noise (incorproates timestep here)

re=rand(NE,1);          ri=rand(NI,1); %connection probabilities?
a=[aE*ones(NE,1);     aI*ones(NI,1)]; %recovery "adaptation current"
b=[bE1*ones(NE/2,1); bE2*ones(NE/2,1);     bI*ones(NI,1)]; %unsure what this is..
tauw=[tauwE*ones(NE,1);     tauwI*ones(NI,1)]; %also for tau..
%C=[CE*ones(NE,0);     CI*ones(NI,1)];
gL=[gLE*ones(NE,1);   gLI*ones(NI,1)]; %leak conductance
%he's stacking exitatiory ontop of inhibitiory cells in this vector

E = ELE; %why are these both for excitatory only??
C = CE;
vT = VTE; %vT is only excitatory here?
vR = [VRE*ones(NE,1);   VRI*ones(NI,1)]; %VT isn't split, but VR is
deltaT = deltaTE;

G = [ WEE0*rand(NE/2,NE/2), zeros(NE/2,NE/2), WIE0*rand(NE/2,NI/2), zeros(NE/2,NI/2); ...
    zeros(NE/2,NE/2), WEE0*rand(NE/2,NE/2), zeros(NE/2,NI/2), WIE0*rand(NE/2,NI/2); ...
    zeros(NI/2,NE/2), WEI0*rand(NI/2,NE/2), WII0*rand(NI/2,NI);
    WEI0*rand(NI/2,NE/2), zeros(NI/2,NE/2),  WII0*rand(NI/2,NI)];
%size is total neurons x total neurons
%wierd part of this matrix is the inhibtory-to-exitatory, scaled much differently (can see with imagesc)
% I think these are connection strengths x random number..
%made up of these We->i etc matricies
%zeros in these matricies are the non We->i connections

%see how G lines up with N cell length vectors, might clarify G

tauSE = 50;
tauSI = 5;

%I = [I0*ones(NE,1); I0*ones(NI,1)];

Ibias_index = 1;
Ibiastrials = 0.1:.02:0.1;
%Ibiastrials = 0.02:0.02:0.02; %test 1 trial only, low Ibias for many switches
Ntrials = length(Ibiastrials);

Nswitches = zeros(1,Ntrials);
switchtime = zeros(1000,Ntrials);
switchofftime = zeros(1000,Ntrials);

T1 = zeros(1,Ntrials);
T2 = zeros(1,Ntrials);

Ibias2 = 0.1; %low Ibias for switches
trial = 0;

%%%short term synaptic plasticity%%% %D1 ONLY%
%p0 = 0.5; %base probability release
tau_d = 9200; %http://bmsr.usc.edu/files/2012/09/JCN1.pdf
tau_d = 5; %http://bmsr.usc.edu/files/2012/09/JCN1.pdf
%might be able to just get rid of tau here, whatever it does..

%tau_f = 0;
%f_f = 0; %facilitation release factor
%f0 = 0;

D = ones(size(tvector));
F = ones(size(tvector));

pr = 0.01; %po*F %1/6
%synaptic depression/gating pr
G(:,1:NE) = G(:,1:NE)/pr; %Only applied to excitatory neurons
E_revE = 0; %reversal potential
E_revI = E;


for Ibias1 = Ibiastrials;
    switchcount = 0;
    trial = trial+1
    
    v = E*ones(NE+NI,length(tvector)); %membrane potential for every cell
    w = zeros(NE+NI,length(tvector)); %connectivity strength? (might be Isra actually)
    v(:,1)=-80+40*rand(NE+NI,1);    % Initial values of v
    w(:,1)= 40*rand(NE+NI,1);        % Initial values of u
    %not sure why this is 40
    
    S = zeros(NE+NI,length(tvector));
    D_cuml = ones(NE+NI, length(tvector)); %avg D over time
    firings = zeros(50*(NE+NI)*tmax/1000,2);
    nspikes = 0;
    nEspikes = 0;
    switchon = 0;
    switchtest = zeros(1,length(tvector));
    switchtestalt = zeros(1,length(tvector));
    switchbias = zeros(1,length(tvector));
    showme = zeros(size(v));
    for j = 1:length(tvector)-1            % simulation of tmax
        
        if ( mod(switchcount,2) == 0 )
            Ibias = Ibias1; %switch between 1 bias and another
        else
            Ibias = Ibias2;
        end
        
        if ( switchon == 1 ) %if switch = true
            %give both pools the same current
            Ibias = 0; % give no bias or switch current
            %Iswitch = 0;
            Iswitch = I0;
        else
            %Iswitch = Iswitch0; %otherwise, give I0
            
            %give both pools the same current
            Ibias = 0; % give no bias or switch current
            Iswitch = I0;
        end
        
        if ( tvector(j) < t_init ) %immediate pulse
            I=[I0 + Isigma*randn(NE/2,1); Iswitch*0.5 + Isigma*randn(NE/2,1);...
                I0I + Isigma*randn(NI,1)]; % thalamic input
        else
            I=[I0*(1+Ibias) + Isigma*randn(NE/2,1); Iswitch + Isigma*randn(NE/2,1);...
                I0I + Isigma*randn(NI,1)]; % thalamic input
        end
        fired=find(v(:,j)>20);% v = 20   % indices of spikes
        %         if ~isempty(fired)
        %             keyboard %gimme some action
        %         end
        Efired = fired(find(fired<NE+1));
        Ifired = fired(find(fired>NE));
        showme(fired,j) = 1;
        Espikes = length(Efired);
        
        spikes = length(fired);
        firings(nspikes+1:nspikes+spikes,:)=[tvector(j)+0*fired,fired];
        
        %update synaptic gating variable --> increase # open channels
        S(Efired,j) = S(Efired, j) + (pr*D_cuml(Efired, j).*(1 - S(Efired, j)));
        %this is how you update synpatic gating for the ones that fired
        
        %update depression after firing, reduction in vesicles
        S(Ifired,j) = 1;
        
        %S(fired,j) = S(fired,j) + 1; %original update (no depression)
        nspikes = nspikes + spikes;
        nEspikes = nEspikes + Espikes;
        v(fired,j)=vR(fired); %after spike
        w(fired,j) = w(fired,j) + b(fired); %after spike
        %not sure why w gets updated by b
        
        I = I + (E_revE - v(:,j)).*(G(:,1:NE)*S(1:NE,j)); %only the excitatory neurons
        %here G dims are 250 x 200
        %total neurons x excitatory neruons
        %matrix seems to be made of these We->i etc matricies
        
        I = I + (E_revI - v(:,j)).*(G(:,NE+1:NE+NI)*S(NE+1:NE+NI,j)); %only the inhibitory neurons
        
        v(:,j+1)=v(:,j)+dt*((E-v(:,j)+deltaT*exp((v(:,j)-vT)/deltaT)).*gL-w(:,j)+I)/C;
        %    w(:,j+1) = w(:,j).*exp(-dt./tauw);
        %         v(:,j+1) = max(v(:,j+1),min(E,min(vR)));
        %         v(:,j+1) = min(v(:,j+1),25);
        
        
        wss = a.*(v(:,j)-E);
        w(:,j+1) = wss + (w(:,j)-wss).*exp(-dt./tauw);
        %    w(:,j+1)=w(:,j)+dt*(a.*(v(:,j+1)-E)-w(:,j))./tauw;
        %update synaptic and depression variables between spikes
        S(:,j+1) = [S(1:NE,j)*exp(-dt/tauSE); S(NE+1:end,j)*exp(-dt/tauSI)];
        
        D_cuml(Efired,j) = D_cuml(Efired,j)*(1-pr);
        D_cuml(:,j+1) = 1 + (D_cuml(:,j) - 1)*exp(-dt/tau_d); %gives a slightly different number than my equation from the book... weird
        switchtestalt(j) = mean(S(1:NE/2,j+1));
        switchtest(j) = mean(S(NE/2+1:NE,j+1));
        %compare excitatory & inhibitory mean S, these are the two groups
        
        
        switchbias(j) = Ibias;
        
        if ( tvector(j) > t_init )
            if  ( (switchtest(j) > 4*switchtestalt(j) ) && (switchon == 0 ) )
                switchcount = switchcount + 1
                switchtime(switchcount,trial) = tvector(j);
                switchon = 1;
            end
            
            if ( (switchon == 1 ) && ( switchtestalt(j) > 4*switchtest(j) ) )
                switchon = 0;
                switchofftime(switchcount,trial) = tvector(j);
                
            end
        end
        
    end;
    Nswitches(trial) = switchcount;
    mean_Erate = nEspikes/(NE*tmax/1000)
    mean_Irate = (nspikes-nEspikes)/(NI*tmax/1000)
    durations = switchtime(2:Nswitches(trial),trial) - ...
        switchofftime(1:Nswitches(trial)-1,trial);
    
    T2s = durations(1:2:end);
    T1s = durations(2:2:end);
    
    trial_save = strcat(duration_file,'_',num2str(Ibias1),'_',num2str(Ibias2),'.mat');
    save(trial_save, 'T1s', 'T2s' );
    
    Nswitches = length(T1s) + length(T2s)
    
    T1(trial) = mean(T1s);
    T2(trial) = mean(T2s);
    nspikes
    format shorte
    data = [Ibias1 Ibias2 T1(trial) T2(trial)]
    
    figure()
    subplot(2,1,1)
    hist(T1s)
    drawnow
    subplot(2,1,2)
    hist(T2s)
    drawnow
end


figure()
plot(Ibiastrials,log(T1))

hold on

plot(Ibiastrials,log(T2),'r')

figure()

plot(Ibiastrials,T2./(T1+T2))

save_file = strcat(basefile,'.mat');
save (save_file,'Ibias*', '*switch*', 'T*', 'Ns*')

close all

spikeplot = make_spikeplot(showme);
figure(1)
imagesc(spikeplot)

poolA = [1:NE/2, NE+1:NE+NI/2];
poolB = [(NE/2)+1:NE, (NE+1)+NI/2:NE+NI];


reordered_plot = NaN(size(spikeplot));
reordered_plot(1:numel(spikeplot(:,1))/2,:) = spikeplot(poolA,:);
reordered_plot((numel(spikeplot(:,1))/2)+1:end,:) = spikeplot(poolB,:);
figure(2)
imagesc(reordered_plot)
xlabel('timepoints')
ylabel('cell')
switch disable_init_pulse
    case 'yes'
        title({'Paul''s model with equal pools & currents','spiking matrix reordered by pool (pool A = rows 1:125)'})
        print('pauls_model_EqualPools_and_Currents','-djpeg')
    case 'no'
        title({'Paul''s model with equal pool & currents with inital pulse','spiking matrix reordered by pool (pool A = rows 1:125)'})
        print('pauls_model_EqualPoolsAndCurrents_and_InitalPulse','-djpeg')
end




