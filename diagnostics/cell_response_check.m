clear
clc
close all
format compact
rng('shuffle') %this is probably important...

%basedir = '~/Desktop/ksander/rotation/project/';
basedir = '~/Desktop/work/ACClab/rotation/project/';
addpath(basedir)
addpath(fullfile(basedir,'helper_functions'))

%----12/19/18: this updates v26 model with slow depression option. 
%Must recalibrate network for:
%   a) fast switching when all things are equal (~1-2 sec)
%   b) stimulus prolongs switching
%   c) spikerates in the I = ~30hz and E = ~15hz sweet zone


tmax = 1.5;
Ispk_time = 1.25; 
Espk_time = .25;
Icell_2_spike = 101;

options = set_options('modeltype','diagnostics','comp_location','bender',...
    'tmax',tmax,'percent_Dslow',.5,...
    'stim_pulse',[tmax,0],'cut_leave_state',tmax,'sample_Estay_offset',0);

options.EtoE = .0405 *1;   %1
options.ItoE = 8; %3
options.EtoI = .35; %2 are best 
options.stim_targs = 'baseline'; %'baseline' | 'Estay' |'baseline'
Rext = 1400 *1; %poisson spike train rate for noise, Hz
Rstim = 300; %rate for stimulus input spikes

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
Erev = 0; %reversal potential, excitatory
Irev = -70e-3; %reversal potential, inhibitory
Gg = 10e-9; %max connductance microSiemens
Pr = NaN(pool_options.num_cells,1); %release probability
Pr(celltype.excit) = .2; %excitatory release probability
Pr(celltype.inhib) = .2; %inhibitory release probability
Td.fast = .3;%synaptic depression time constant, seconds (fast)
Td.slow = 7;%slow time constant
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
Tau_ext(celltype.excit) = 3.5e-3; % was 2e-3
Tau_ext(celltype.inhib) = 2e-3; % was 5e-3; 
initGext = 10e-9; %noisy conductance initialization value, nano Siemens
deltaGext = 1e-9; %increase noisy conducrance, nano Siemens
%Rext = 1400; %poisson spike train rate for noise, Hz
%---adaptation conductance----
Vth = -50e-3; %ALEIF spike threshold mV
delta_th = 2e-3; %max voltage threshold, mV  (deltaVth in equation)
Tsra = NaN(pool_options.num_cells,1);%adaptation conductance time constant, ms
Tsra(celltype.excit) = 25e-3; %excitatory
Tsra(celltype.inhib) = 25e-3; %inhibitory
detlaGsra = 12.5e-9; %increase adaptation conductance, nano Siemens
%----timecourse---------------
timestep = options.timestep;
timevec = 0:timestep:options.tmax;
num_timepoints = numel(timevec);
Lext = Rext * timestep; %poisson lambda for noisy conductance

%preallocate variables
%-------------------------------------------------------------------------
%---stimuli info--------------
stim_info = struct();
switch options.stim_targs
    case 'Eswitch'
        stim_info.targ_cells = celltype.excit & celltype.pool_switch; %Eswitch cells
    case 'Estay'
        stim_info.targ_cells = celltype.pool_stay & celltype.excit; %Estay
    case 'baseline'
        stim_info.targ_cells = logical(zeros(pool_options.num_cells,1)); %no targets
end
stim_info.stimA_lambda = Rstim * timestep; %poisson lambda for stimulus conductance
stim_info.stimB_lambda = Rstim * timestep;
stim_info.num_cells = pool_options.num_cells; %just so I don't have to pass pool_options as well
if all(~isnan(options.stim_pulse))
    stim_info.delivery = 'pulse';
    stim_info.pulse = options.stim_pulse ./ timestep;
    stim_info.sample_schedule = options.stim_schedule;
else
    stim_info.delivery = 'constant';
end
%---membrane potential--------
V = NaN(pool_options.num_cells,2);
V(:,1) = El; %inital value of membrane potential is leak potential
%---noisy conductance---------
Gext = NaN([size(V),2]); %noisy conductance (do I & E input in 3rd D)
ext_inds.I = 1;
ext_inds.E = 2;
Gext(:,1,:) = initGext; %initialize at leak conductance
%---adaptation conductance----
Gsra = NaN(size(V));
Gsra(:,1) = 0;
%---gating & depression-------
Sg = NaN(size(V)); %synaptic gating
Sg(:,1) = 0; %initalize at zero??
Dmax.fast = 1 - options.percent_Dslow;
Dmax.slow = options.percent_Dslow;
Dfast = NaN(size(V)); %synaptic depression: fast 
Dslow = NaN(size(V)); %synaptic depression: slow 
Dfast(:,1) = Dmax.fast; %initalize at ratio of slow/fast vessicles
Dslow(:,1) = Dmax.slow;
%---spikes--------------------
spikes = zeros(pool_options.num_cells,num_timepoints); %preallocating the whole thing in this one...
%---state tracker-------------
durations = {}; %record duration time, state/stimulus label
state = init_statevar(celltype,options);
options.sample_Estay_offset = options.sample_Estay_offset / timestep;
state.now = state.undecided; %this will always be true when V init to El
%---last init-----------------
experiment_set2go = false; %skip bistability check 
avail_stim = true; %in between stimulus delivery pulses in stay state
timepoint_counter = 1;
idx = 2; %keep indexing vars with idx fixed at 2

fprintf(':::Starting simulation:::\n')

while timepoint_counter < num_timepoints
    
    timepoint_counter = timepoint_counter+1;
    state.timeidx = timepoint_counter; %just so I don't have to pass a million things...
    
    %loop equations
    
    %%excitatiory input
    %I = (Erev - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*Gg;
    %%inhibitory input 
    %I = I + (Irev - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*Gg;
    inp_Ecells = (Erev - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*Gg;
    inp_Icells = (Irev - V(:,idx-1)).*(W(:,celltype.inhib)*Sg(celltype.inhib,idx-1)).*Gg;
    I = inp_Ecells + inp_Icells;
    I = I + (Gext(:,idx-1,ext_inds.E).*(Erev-V(:,idx-1)));
    I = I + (Gext(:,idx-1,ext_inds.I).*(Irev-V(:,idx-1)));
    dVdt = ((El-V(:,idx-1)+(delta_th.*exp((V(:,idx-1)-Vth)./delta_th)))./Rm)...
        + (Gsra(:,idx-1).*(Ek-V(:,idx-1))) + I;
    
    V(:,idx) = ((dVdt./Cm) .* timestep) + V(:,idx-1);
        
    tGext = squeeze(Gext(:,idx-1,:)); %flattened t-1 Gext 
    Gext(:,idx,:) = tGext - ((tGext./Tau_ext) .* timestep); %noisy conductance
    Gsra(:,idx) = Gsra(:,idx-1) - ((Gsra(:,idx-1)./Tsra) .* timestep); %adaptation conductance
    Dfast(:,idx) = Dfast(:,idx-1) + (((Dmax.fast - Dfast(:,idx-1))./Td.fast) .* timestep); %fast syn. depression
    Dslow(:,idx) = Dslow(:,idx-1) + (((Dmax.slow - Dslow(:,idx-1))./Td.slow) .* timestep); %slow syn. depression
    Sg(:,idx) = Sg(:,idx-1) - ((Sg(:,idx-1)./Tsyn) .* timestep); %synaptic gating
    
    spiking_cells = V(:,idx) > spike_thresh;
    
    if timepoint_counter == Espk_time / timestep
        spiking_cells(10) = 1; %single excitatory spike
    elseif timepoint_counter == Ispk_time / timestep
        spiking_cells(Icell_2_spike) = 1; %single inhibitory spike 
    end
    %if sum(spiking_cells) > 0
    %    spiking_cells = V(:,idx) > spike_thresh;
        Gsra(spiking_cells,idx) = Gsra(spiking_cells,idx) + detlaGsra; %adaptation conductance
        Pr_spike = Pr(spiking_cells); 
        %vessicle release for slow/fast vessicles
        fast_release = Pr_spike .* Dfast(spiking_cells,idx); 
        slow_release = Pr_spike .* Dslow(spiking_cells,idx); 
        %synaptic gating, depends on combined vessicle release 
        Sg(spiking_cells,idx) = Sg(spiking_cells,idx) + ...
            ((fast_release + slow_release).*(1-Sg(spiking_cells,idx))); 
        %depression update
        Dfast(spiking_cells,idx) = Dfast(spiking_cells,idx) - fast_release; 
        Dslow(spiking_cells,idx) = Dslow(spiking_cells,idx) - slow_release; 
        V(spiking_cells,idx) = Vreset;
        spikes(spiking_cells,timepoint_counter) = 1;
    %end
    
%     %test for state transition & determine stim availability
%     if experiment_set2go
%         [state,durations] = test4switch(Sg(:,idx),state,durations);
%         [state,avail_stim] = check_stim_avail(stim_info,state);
%     else
%         state.count = state.count + 1; %if you don't run test4switch(), must update this counter outside
%     end
    
%     %run the bistability check
%     if ~experiment_set2go %during bistability check, check_bistability() handles pulse input spikes
%         [BScheck,Pspikes,state] = check_bistability(Sg(:,idx),state);
%         %add pulse spikes (same as below)
%         Gext(:,idx,ext_inds.E) = Gext(:,idx,ext_inds.E) + (deltaGext.*Pspikes);
%         switch BScheck.status
%             case 'fail'
%                 fprintf(':::Bistability check failure:::')
%                 fprintf('---at t=%.2f(s)',timepoint_counter*timestep)
%                 return
%             case 'pass'
%                 fprintf('---passed bistability check')
%                 experiment_set2go = true; %we're ready to roll
%         end
%     end
        
    %input spikes: noise & stimulus
    %---noisy spiking input from elsewhere
    %ext_spikes = poissrnd(Lext,pool_options.num_cells,2);
    ext_spikes = zeros(pool_options.num_cells,2);
    if timepoint_counter > Ispk_time / timestep - (100e-3 / timestep)
        ext_spikes(celltype.excit,:) = poissrnd(Lext,sum(celltype.excit),2);
    end

%     if ~avail_stim
%         %don't do this for current testing
%         %         %we're in between stimulus delivery pulses
%         %         %do half-noise to all E-cells-- gives barely spiking undecided state
%         %         ext_spikes(celltype.excit,:) = poissrnd(Lext*.5,sum(celltype.excit),2); %half noise E-cells
%     elseif avail_stim && strcmp(stim_info.delivery,'pulse') && experiment_set2go
%         %stimulus is available, and we're doing pulse-sample delivery
%         Tsample = sum(stim_info.pulse); %how long for a single on, off sequence
%         Tsample = mod(state.sample_clock,Tsample); %find out how far into the sample
%         if Tsample <= options.sample_Estay_offset
%             %if during first Xms of a sample, give E-stay full noise &
%             %E-switch half-noise to kick on the stay-state
%             ext_spikes(celltype.excit & celltype.pool_switch,:) = ...
%                 poissrnd(Lext*.5,sum(celltype.excit & celltype.pool_switch),2);
%         end
%     end
%     %---spiking input from stimulus
%     if experiment_set2go
%         stim_spikes = timepoint_stimulus(stim_info,state); %get stimulus spikes
%         %always exitatory, add 'em both together for one calculation
%         ext_spikes(:,ext_inds.E) = ext_spikes(:,ext_inds.E) + stim_spikes;
%     end
    

    %update Gexternal. Don't have to index, they get an increase or zero
    Gext(:,idx,:) = squeeze(Gext(:,idx,:)) + (deltaGext.*ext_spikes);
    
    
    %for recording vars
    Vrec(:,timepoint_counter) = V(:,idx);
    Drec_fast(:,timepoint_counter) = Dfast(:,idx);
    Drec_slow(:,timepoint_counter) = Dslow(:,idx);
    Srec(:,timepoint_counter) = Sg(:,idx);
    Irec(:,timepoint_counter) = inp_Icells;
    
    %lag equation vars for next timepoint
    V = next_timepoint(V);
    Gsra = next_timepoint(Gsra);
    Gext = next_timepoint(Gext);
    Dfast = next_timepoint(Dfast);
    Dslow = next_timepoint(Dslow);
    Sg = next_timepoint(Sg);
    
    
    %check timeout for non-switching
    if state.count >= state.noswitch_timeout
        fprintf(':::Bistability check failure:::')
        TOF = timepoint_counter*timestep;
        fprintf('---no switch timeout at t=%.2f(s)',TOF)
        return
    end
    
    
    %progress tracking...
    if mod(timepoint_counter,floor(num_timepoints * .10)) == 0 %10 percent
        progress = (timepoint_counter /  num_timepoints) * 100;
        fprintf('Simulation %.1f percent complete\n',progress);
    end
    
end



lnsz = 3; %spikerate plots
fontsz = 12;




figure()
imagesc(Vrec)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);

Yticks = num2cell(get(gca,'YLim'));
Yticks = Yticks{2};
Yticks = [Yticks*.25,Yticks*.75];
Ylabs = {'stay pool','leave pool'};
ytickangle(90)
set(gca,'Xdir','normal','Ytick',Yticks,'YTickLabel', Ylabs);
xlabel('time (s)','FontWeight','b')
title('Vm')
colb = colorbar;


Vresps = max(Vrec(:,Espk_time/ timestep:end),[],2);
%look for inhibitiory cells 
Vmax = max(Vresps(celltype.inhib));
Vmax = find(Vresps == Vmax & celltype.inhib,1);

figure()
plot(timevec ,Vrec(Vmax,:) ./ 1e-3,'LineWidth',2)
title('membrane voltage')
xlabel('time (s)')
ylabel('mV')


inib_inp = min(Irec(:,Ispk_time / timestep:end),[],2);
%look for excitatory cells 
Imin = min(inib_inp(celltype.excit));
Imin = find(inib_inp == Imin & celltype.excit,1);

figure()
Tminus = (Ispk_time / timestep) - (10e-3 / timestep); %10 ms before timestep
Tend = (Ispk_time / timestep) + (100e-3 / timestep); %100 ms after spike 
plot(timevec(Tminus:Tend) ./ 1e-3,Irec(Imin,Tminus:Tend) ./ 1e-12,'LineWidth',2)
title('inhibitory input on E-cells')
xlabel('T+ spike (ms)')
ylabel('current (pico amps)')
axis tight

figure()
%do the raster plot
spikeplot = make_spikeplot(spikes);
%reorganize the cell groups
raster = [spikeplot(celltype.pool_stay & celltype.excit,:);...
    spikeplot(celltype.pool_switch & celltype.inhib,:);...
    spikeplot(celltype.pool_switch & celltype.excit,:);...
    spikeplot(celltype.pool_stay & celltype.inhib,:)];
imagesc(raster)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);

Yticks = num2cell(get(gca,'YLim'));
Yticks = Yticks{2};
Yticks = [Yticks*.25,Yticks*.75];
Ylabs = {'stay pool','leave pool'};
ytickangle(90)
set(gca,'Xdir','normal','Ytick',Yticks,'YTickLabel', Ylabs);
xlabel('time (s)','FontWeight','b')
title({'spikes','(spikes in matrix enlarged for visualization)'})

figure()
zoomspikes = zeros(size(spikes));
zoomspikes(Icell_2_spike,Ispk_time / timestep) = 1;
zoomspikes = zoomspikes + spikes;
raster = [zoomspikes(celltype.pool_stay & celltype.excit,:);...
    zoomspikes(celltype.pool_switch & celltype.inhib,:);...
    zoomspikes(celltype.pool_switch & celltype.excit,:);...
    zoomspikes(celltype.pool_stay & celltype.inhib,:)];
imagesc(raster(:,Tminus:Tend))
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',...
    ((x - (Ispk_time / timestep - Tminus))*timestep)/1e-3),...
    Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
Yticks = num2cell(get(gca,'YLim'));
Yticks = Yticks{2};
Yticks = [Yticks*.25,Yticks*.75];
Ylabs = {'stay pool','leave pool'};
ytickangle(90)
set(gca,'Xdir','normal','Ytick',Yticks,'YTickLabel', Ylabs);
xlabel('T+ spike (ms)','FontWeight','b')
%title({'spikes','(spikes in matrix enlarged for visualization)'})



window_sz = 50e-3;
Dmu_fast = sim_windowrate(Drec_fast,timestep,celltype,window_sz);
Dmu_fast = structfun(@(x) x.* timestep ,Dmu_fast,'UniformOutput',false); %undo hz conversion

Dmu_slow = sim_windowrate(Drec_slow,timestep,celltype,window_sz);
Dmu_slow = structfun(@(x) x.* timestep ,Dmu_slow,'UniformOutput',false); %undo hz conversion

figure;
%plot the aggregated timecourses
subplot(2,1,1);hold on
plot(Dmu_fast.Estay,'Linewidth',lnsz)
plot(Dmu_fast.Eswitch,'Linewidth',lnsz)
plot(Dmu_fast.Istay,'Linewidth',lnsz)
plot(Dmu_fast.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Fast depression')
ylabel({sprintf('Pool average  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz);axis tight;hold off
subplot(2,1,2);hold on
plot(Dmu_slow.Estay,'Linewidth',lnsz)
plot(Dmu_slow.Eswitch,'Linewidth',lnsz)
plot(Dmu_slow.Istay,'Linewidth',lnsz)
plot(Dmu_slow.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('slow depression')
ylabel({sprintf('Pool average  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
set(gca,'FontSize',fontsz);axis tight;hold off

%spikerates by sliding window...
Srate = sim_windowrate(spikes,timestep,celltype,window_sz);

figure;
%plot the aggregated timecourses
hold on
plot(Srate.Estay,'Linewidth',lnsz)
plot(Srate.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight;hold off
keyboard

figure;
%plot the aggregated timecourses
hold on
plot(Srate.Estay,'Linewidth',lnsz)
plot(Srate.Eswitch,'Linewidth',lnsz)
plot(Srate.Istay,'Linewidth',lnsz)
plot(Srate.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight;hold off


figure;
plot(Srate.Estay-Srate.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Spiking')
ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight

Sg_mu.Estay = mean(Srec(celltype.excit & celltype.pool_stay,:),1);
Sg_mu.Eswitch = mean(Srec(celltype.excit & celltype.pool_switch,:),1);
Sg_mu.Istay = mean(Srec(celltype.inhib & celltype.pool_stay,:),1);
Sg_mu.Iswitch = mean(Srec(celltype.inhib & celltype.pool_switch,:),1);

%plot the aggregated timecourses
figure;
hold on
plot(Sg_mu.Estay,'Linewidth',lnsz)
plot(Sg_mu.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
%legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
legend({'E-stay','E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight;hold off


%plot the aggregated timecourses
figure;
hold on
plot(Sg_mu.Estay,'Linewidth',lnsz)
plot(Sg_mu.Eswitch,'Linewidth',lnsz)
plot(Sg_mu.Istay,'Linewidth',lnsz)
plot(Sg_mu.Iswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight;hold off


%plot the aggregated timecourses
figure;
plot(Sg_mu.Estay-Sg_mu.Eswitch,'Linewidth',lnsz);hold on
plot(Sg_mu.Istay-Sg_mu.Iswitch,'Linewidth',lnsz);hold off
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch','I-stay minus I-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight

figure
plot(Sg_mu.Estay-Sg_mu.Eswitch,'Linewidth',lnsz)
Xticks = num2cell(get(gca,'Xtick'));
Xlabs = cellfun(@(x) sprintf('%.1f',x*timestep),Xticks,'UniformOutput', false); %this is for normal stuff
set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
title('Synaptic gating')
ylabel({sprintf('Mean pool S-gating  (%ims bins)',window_sz*1e3)})
xlabel('time (s)')
legend({'E-stay minus E-switch'},'location','northoutside','Orientation','horizontal')
set(gca,'FontSize',fontsz)
axis tight

if numel(durations) > 1
    timecourse = size(durations);
    timecourse(2) = timecourse(2) + 1;
    timecourse = cell(timecourse);
    timecourse(:,1:3) = cellfun(@(x) x*options.timestep,durations(:,1:3),'UniformOutput',false);
    timecourse(:,3) = cellfun(@(x) mod(x,sum(options.stim_pulse)),timecourse(:,3),'UniformOutput',false);
    %current sample's onset, rounding is needed for subsequent operations
    timecourse(:,4) = cellfun(@(x,y) round(x-y,2),timecourse(:,1),timecourse(:,3),'UniformOutput',false);
    samp_onsets = unique(cat(1,timecourse{:,4})); %like unique won't work properly here without rounding
    timecourse(:,4) = cellfun(@(x) find(x==samp_onsets),timecourse(:,4),'UniformOutput',false); %would also break without rounding
    timecourse(:,end) = durations(:,end);
    %timecourse(:,1:end-1) = cellfun(@(x) sprintf('%.3f',x),timecourse(:,1:end-1),'UniformOutput',false); %for printing
    timecourse = cell2table(timecourse,'VariableNames',{'event_time','duration','sample_time','sample_number','state'});
    figure
    state_durs = timecourse.duration(startsWith(timecourse.state,'stim'));
    if numel(state_durs) > 0
        if numel(state_durs) < 15
            histogram(state_durs,numel(state_durs))
        else
            histogram(state_durs)
        end
        title(sprintf('stay-state durations\nmu = %.2f',mean(state_durs)))
        ylabel('seconds');xlabel('frequecy');set(gca,'FontSize',fontsz)
    end
end


function cell_data = sim_spikerate(cell_raster,timestep,celltype)

num_binsamps = 75e-3/timestep; %num samples in 2ms
raster_sz = size(cell_raster);
if mod(raster_sz(2),num_binsamps) ~= 0 %you have to trim it down, equally divisible by bin size
    cell_raster = cell_raster(:,1:end - mod(raster_sz(2),num_binsamps));
    raster_sz = size(cell_raster); %should be good now
end
bin_magic = [raster_sz(1), num_binsamps,raster_sz(2)/num_binsamps]; %set up for a magic trick
cell_raster = reshape(cell_raster,bin_magic);
cell_raster = squeeze(sum(cell_raster,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_raster = repmat(cell_raster,[num_binsamps 1 1]);
cell_raster = reshape(cell_raster,raster_sz); %put the rabit back in the hat

%normal people indexing that makes sense, then take the mean
cell_data.Estay = mean(cell_raster(celltype.excit & celltype.pool_stay,:),1);
cell_data.Eswitch = mean(cell_raster(celltype.excit & celltype.pool_switch,:),1);
cell_data.Istay = mean(cell_raster(celltype.inhib & celltype.pool_stay,:),1);
cell_data.Iswitch = mean(cell_raster(celltype.inhib & celltype.pool_switch,:),1);
end

function cell_data = sim_windowrate(cell_raster,timestep,celltype,window_sz)

num_binsamps = window_sz/timestep; %num samples in Xms
cell_raster = num2cell(cell_raster,2);
k = ones(1, num_binsamps);
k = k ./ (num_binsamps * timestep); %convert to Hz
cell_raster = cellfun(@(x) conv(x, k, 'same'),cell_raster,'UniformOutput',false);
cell_raster = cat(1,cell_raster{:});

%normal people indexing that makes sense, then take the mean
cell_data.Estay = mean(cell_raster(celltype.excit & celltype.pool_stay,:),1);
cell_data.Eswitch = mean(cell_raster(celltype.excit & celltype.pool_switch,:),1);
cell_data.Istay = mean(cell_raster(celltype.inhib & celltype.pool_stay,:),1);
cell_data.Iswitch = mean(cell_raster(celltype.inhib & celltype.pool_switch,:),1);
end



