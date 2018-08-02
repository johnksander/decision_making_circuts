%LIF_model_sra.m
clear

%based off chapter 1 tutorial 1c: Models based on extensions of the LIF neuron.

figure(1)
clf

tau = 0.010;
dt = 0.0001;
t = 0:dt:0.5;
ton = 0.15;
toff = 0.35;
non = round(ton/dt);
noff = round(toff/dt);
tref = 0.002;
tref_method = 1;

E_L = -0.070;
E_K = -0.080;

Vth0 = -0.050;
Vreset = -0.080;
Cm = 100e-12;
G_L = Cm/tau;

tsra = 0.200;
delta_G = 1e-9;
Iapp = [240e-12 400e-12];
Ntrials = length(Iapp);

for trial = 1:Ntrials;
    I = zeros(size(t));
    I(non:noff) = Iapp(trial)*ones(1,noff+1-non);
    V = E_L*ones(size(t));
    spikes = zeros(size(t));
    Gsra = zeros(size(t));
    
    Vth = Vth0*ones(size(t));
    t_last_spike = -10*tref;
    
    for i = 2:length(t)
        Gsra(i) = Gsra(i-1)*exp(-dt/tsra);
        Vth(i) = Vth0 + (Vth(i-1)-Vth0)*exp(-dt/tref);
        Vss = ( I(i) + G_L*E_L + Gsra(i)*E_K)/(G_L + Gsra(i));
        taueff = Cm/(G_L+Gsra(i));
        V(i) = Vss + ( V(i-1)-Vss)*exp(-dt/taueff);
        if ( tref_method == 1 )
            if ( t(i) < t_last_spike + tref )
                V(i) = Vreset;
            end
        end
        
        if V(i) > Vth(i)
            spikes(i) = 1;
            Gsra(i) = Gsra(i) + delta_G;
            V(i) = Vreset;
        end
    end
    
    figure(1)
    hold on
    subplot('position',[0.44*trial-0.27 0.56 0.36 0.36])
    plot(t,V,'k');
    if ( trial == 1 )
        ylabel('V_m (V)')
        set(gca,'YTick',[-0.08 -0.04 0])
    else
        set(gca,'YTickLabel','')
    end
    axis([0 0.5 Vreset-0.005 Vth0+0.005])
    set(gca,'XTickLabel','')
    if ( trial == 1)
        title('240pA')
    end
    if ( trial == 2)
        title('400pA')
    end
    subplot('position',[0.44*trial-0.27 0.13 0.36 0.36])
    plot(t,Gsra*1e9,'k')
    xlabel('Time (sec)')
    if ( trial == 1 )
        ylabel('G_{SRA} (n\Omega)')
    else
        set(gca,'YTickLabel','')
    end
    
    axis([0 0.5 0 8e9*delta_G])
    
    %         set(gca,'XTick',[0.26 0.28])
    %         set(gca,'XTickLabel',{'0.26', '0.28'})
    
    
end
