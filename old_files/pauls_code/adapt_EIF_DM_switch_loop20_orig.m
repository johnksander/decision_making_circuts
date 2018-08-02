% adapting_DM.m
% decision-making network with adaptation to implement switching
clear

basefile = 'DM_switch_STP_20';
    
duration_file = strcat(basefile, 'dur');
    

set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
NE=200; %excitatory neurons
NI=50; %inhibitory neurons
WEE0 = 5/NE;
WEI0 = 20/NE; %stronger e- i than e-e
WIE0 = 300/NI;%inhibitory - excitatory strongest %*6
WII0 = 0/NI; %inhibitory - inhibitory - no connections
% WEE0 = 0/NE;
% WEI0 = 0/NE;
% WIE0 = 0/NI;
% WII0 = 0/NI;

aE= 2; %time scale of recovery (typical = .02, so this is really fast)
aI = 2; %time scale of recovery (very fast)
CE = 200;
CI = 200;
gLE = 12;
gLI = 10;
ELE = -62;
ELI = -62;
VTE = -50; 
VTI = -50;
deltaTE = 2;
deltaTI = 2;
tauwE = 300;
tauwI = 30;
bE1 = 3;
bE2 = 0.5;
bI = 0;
VRE = -58;
VRI = -56;

tmax = 900000;
%tmax = 300000;
dt = .25;   % in ms
tvector = dt:dt:tmax;
t_init = 200;

Nwindow = floor(10000/dt);

I0 = 100;
I0I = 80;
Iswitch0 = 120;

Isigma = 400/sqrt(dt);

re=rand(NE,1);          ri=rand(NI,1);
a=[aE*ones(NE,1);     aI*ones(NI,1)];
b=[bE1*ones(NE/2,1); bE2*ones(NE/2,1);     bI*ones(NI,1)];
tauw=[tauwE*ones(NE,1);     tauwI*ones(NI,1)];
%C=[CE*ones(NE,0);     CI*ones(NI,1)];
gL=[gLE*ones(NE,1);   gLI*ones(NI,1)];
E = ELE;
C = CE;
vT = VTE;
vR = [VRE*ones(NE,1);   VRI*ones(NI,1)];
deltaT = deltaTE;

G = [ WEE0*rand(NE/2,NE/2), zeros(NE/2,NE/2), WIE0*rand(NE/2,NI/2), zeros(NE/2,NI/2); ...
    zeros(NE/2,NE/2), WEE0*rand(NE/2,NE/2), zeros(NE/2,NI/2), WIE0*rand(NE/2,NI/2); ...
    zeros(NI/2,NE/2), WEI0*rand(NI/2,NE/2), WII0*rand(NI/2,NI);
    WEI0*rand(NI/2,NE/2), zeros(NI/2,NE/2),  WII0*rand(NI/2,NI)];

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

%tau_f = 0;
%f_f = 0; %facilitation release factor
%f0 = 0;

D = ones(size(tvector));
F = ones(size(tvector));

pr = 0.01; %po*F %1/6
G(:,1:NE) = G(:,1:NE)/pr;
E_revE = 0;
E_revI = E;


for Ibias1 = Ibiastrials;
    switchcount = 0;
    trial = trial+1
    
    v = E*ones(NE+NI,length(tvector));
    w = zeros(NE+NI,length(tvector));
    v(:,1)=-80+40*rand(NE+NI,1);    % Initial values of v
    w(:,1)= 40*rand(NE+NI,1);                 % Initial values of u
    
    S = zeros(NE+NI,length(tvector));
    D_cuml = ones(NE+NI, length(tvector)); %avg D over time
    firings = zeros(50*(NE+NI)*tmax/1000,2);
    nspikes = 0;
    nEspikes = 0;
    switchon = 0;
    switchtest = zeros(1,length(tvector));
    switchtestalt = zeros(1,length(tvector));
    switchbias = zeros(1,length(tvector));
    
    for j = 1:length(tvector)-1            % simulation of tmax
        
        if ( mod(switchcount,2) == 0 )
            Ibias = Ibias1; %switch between 1 bias and another
        else
            Ibias = Ibias2;
        end
        
        if ( switchon == 1 ) %if switch = true
            Ibias = 0; % give no bias or switch current
            Iswitch = 0;
        else
            Iswitch = Iswitch0; %otherwise, give I0
        end
        
        if ( tvector(j) < t_init ) %immediate pulse
            I=[I0 + Isigma*randn(NE/2,1); Iswitch*0.5 + Isigma*randn(NE/2,1);...
                I0I + Isigma*randn(NI,1)]; % thalamic input
        else
            I=[I0*(1+Ibias) + Isigma*randn(NE/2,1); Iswitch + Isigma*randn(NE/2,1);...
                I0I + Isigma*randn(NI,1)]; % thalamic input
        end
        fired=find(v(:,j)>20);% v = 20   % indices of spikes
        Efired = fired(find(fired<NE+1));
        Ifired = fired(find(fired>NE));
        
        Espikes = length(Efired);
        
        spikes = length(fired);
        firings(nspikes+1:nspikes+spikes,:)=[tvector(j)+0*fired,fired];
        
        %update synaptic gating variable --> increase # open channels
        S(Efired,j) = S(Efired, j) + (pr*D_cuml(Efired, j).*(1 - S(Efired, j)));
        %update depression after firing, reduction in vesicles
        S(Ifired,j) = 1;
        
        %S(fired,j) = S(fired,j) + 1; %original update (no depression)
        nspikes = nspikes + spikes;
        nEspikes = nEspikes + Espikes;
        v(fired,j)=vR(fired); %after spike
        w(fired,j) = w(fired,j) + b(fired); %after spike

        I = I + (E_revE - v(:,j)).*(G(:,1:NE)*S(1:NE,j));
        I = I + (E_revI - v(:,j)).*(G(:,NE+1:NE+NI)*S(NE+1:NE+NI,j));

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
        D_cuml(:,j+1) = 1 + (D_cuml(:,j) - 1)*exp(-dt/tau_d);         
        switchtestalt(j) = mean(S(1:NE/2,j+1));
        switchtest(j) = mean(S(NE/2+1:NE,j+1));
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


