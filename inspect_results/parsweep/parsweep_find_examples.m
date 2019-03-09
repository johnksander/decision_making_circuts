clear
clc
format compact
close all

%this is broken somehow...

%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_fastD_Rlim_baseline';
basedir = '/home/acclab/Desktop/ksander/rotation/project';
figdir = fullfile(basedir,'Results',['figures_' sim_name]);
resdir = fullfile(basedir,'Results',sim_name);
addpath(fullfile(basedir,'helper_functions'))

%get results file
result_data = load(fullfile(resdir,'summary_file'));
result_data = result_data.result_data;

%get data-of-interest
ItoE = cellfun(@(x)  x.ItoE,result_data(:,2));
EtoI = cellfun(@(x)  x.EtoI,result_data(:,2));
%state_dur = cellfun(@(x) mean(x),result_data(:,1)); %"durations" reserved for output file 
state_dur = cellfun(@(x) mean(log10(x+eps)),result_data(:,1));
Sdurs = 10.^state_dur; %in seconds 

%slow_range = [9.5, 10.5]; %seconds
slow_range = log10([20, 30]); %seconds
fast_range = log10([1, 1.5]);
%many more slow, than fast. So start with the slow networks
in_range = @(x,y) x >= y(1) & x <= y(2);     %find durations X within target range Y (two element vec)


slow_nets = in_range(state_dur,slow_range);
fast_nets = in_range(state_dur,fast_range);

fprintf('\n::::::::slow networks::::::::\n')
fprintf('---range = %.2f - %.2f\n',slow_range)
fprintf('---total networks = %i\n',sum(slow_nets))

fprintf('\n::::::::fast networks::::::::\n')
fprintf('---range = %.2f - %.2f\n',fast_range)
fprintf('---total networks = %i\n',sum(fast_nets))



HM = openfig(fullfile(figdir,'logmean_duration','heatmap.fig'));
hold on
HMdata = gca;

HMdims = get(HMdata,'Children');
HMdims = size(HMdims.CData);
Nx = HMdims(1); Ny = HMdims(2);
Yax = linspace(max(EtoI),min(EtoI),Ny);
Xax = linspace(min(ItoE),max(ItoE),Nx);

nearest_ind = @(x,y) find(abs(x-y) == min(abs(x-y)),1); %return index of closest value 
nearest_val = @(x,y) y(nearest_ind(x,y)); %return closest value 

netX = num2cell(ItoE);
netX = cellfun(@(x) nearest_ind(x,Xax),netX);

netY = num2cell(EtoI);
netY = cellfun(@(x) nearest_ind(x,Yax),netY);

cand_nets = slow_nets | fast_nets;

scatter(netX(cand_nets),netY(cand_nets),'red')
hold off


pause(1)
SP = openfig(fullfile(figdir,'logmean_duration','surface_plot.fig'));pause(1);gcf
hold on

scatter3(ItoE(cand_nets),EtoI(cand_nets),state_dur(cand_nets),40,'black','filled',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1); %mark 'em 

figure()
subplot(2,2,[1:2])
histogram(Sdurs(slow_nets))
title('slow network candidates');xlabel('durations (s)')
subplot(2,2,3)
scatter(EtoI(slow_nets),Sdurs(slow_nets));ylabel('durations (s)');xlabel('EtoI')
subplot(2,2,4)
scatter(ItoE(slow_nets),Sdurs(slow_nets));ylabel('durations (s)');xlabel('ItoE')
keyboard

%test inds 
param_perc = @(x,y) (range(y)*x )+ min(y); %return the Xth percentile value from parameter range Y  
param_perc(.5,ItoE)



net_pair = array2table(NaN(2,3),'VariableNames',{'ItoE','EtoI','duration'},'RowNames',{'slow','fast'});

curr_net = param_perc(1,ItoE(slow_nets));
curr_net = nearest_ind(curr_net,ItoE);

net_pair{'slow','ItoE'} = ItoE(curr_net);
net_pair{'slow','EtoI'} = EtoI(curr_net);
net_pair{'slow','duration'} = Sdurs(curr_net);

%find the matching fast one 
curr_net = nearest_val(net_pair{'slow','ItoE'},ItoE(fast_nets));
curr_net = find(ItoE == curr_net & fast_nets);

net_pair{'fast','ItoE'} = ItoE(curr_net);
net_pair{'fast','EtoI'} = EtoI(curr_net);
net_pair{'fast','duration'} = Sdurs(curr_net);






% scatter3(ItoE(cand_nets),EtoI(cand_nets),log10(state_dur(cand_nets)),20,'black','filled',...
%     'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);pause(1) %mark 'em 


%plot3(ItoE,EtoI,outcome,'.','MarkerSize',15,'color','blue') %give them outlines
%plot3(ItoE,EtoI,outcome,'.','MarkerSize',10,'color','green') 

% 
% scatter3(Nplots(:,1),Nplots(:,2),Nplots(:,3),1000,'red','p','filled',...
%     'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1); %mark 'em w/ stars



%network pair 1:
%slow------ 75 hz
% ItoE      EtoI   duration
%1.2904    0.1948  143.4682
%fast------ 100 hz
% ItoE      EtoI   duration
%1.2877    0.1734    2.1808


%netowork pair 2:
%slow------ 49.7632 hz
% ItoE      EtoI   duration
%0.8069    0.2067  187.2613
%fast------ 59.2903 hz
% ItoE      EtoI   duration
%0.7936    0.1878    2.0959


%network pair 3:
%slow------ 30.0436 hz
% ItoE      EtoI   duration
%0.3679    0.2737  189.1735
%fast------ 175 hz
% ItoE      EtoI   duration
%0.3161    0.2482    2.0199

%network pair 4:
%slow------ 47.0955
% ItoE      EtoI   duration
%0.2800    0.4228  165.4364
%fast------ 175
% ItoE      EtoI   duration
%0.2355    0.4250    2.1626


%network pair 5:
%slow------ 45 hz
% ItoE      EtoI   duration
%0.2921    0.6927  179.6920
%fast------ 175 Hz
% ItoE      EtoI   duration
%0.2119    0.6799    2.2800

%make this into a data structure...

%original --------------------
network_pairs = cell(5,1);
%...
network_pairs{1} = [{1.2904, 0.1948, 75, 'Eswitch'};
    {1.2877,  0.1734, 100, 'Estay'}];

network_pairs{2} = [{0.8069, 0.2067, 49.7632, 'Eswitch'};
    {0.7936, 0.1878, 59.2903, 'Estay'}];

network_pairs{3} = [{0.3679, 0.2737, 30.0436, 'Eswitch'};
    {0.3161, 0.2482, 175, 'Estay'}];

network_pairs{4} = [{0.2800, 0.4228, 47.0955, 'Eswitch'};
    {0.2355, 0.4250, 175, 'Estay'}];

network_pairs{5} =  [{0.2921, 0.6927, 45, 'Eswitch'};
    {0.2119, 0.6799, 175, 'Estay'}];
%end originial --------------------

%now trying to equate------------------
%03212018: I'm just halfing each E-switch stim 
network_pairs{1} = [{1.2904, 0.1948, 75/2, 'Eswitch'};
    {1.2877,  0.1734, 100, 'Estay'}];

network_pairs{2} = [{0.8069, 0.2067, 49.7632/2, 'Eswitch'};
    {0.7936, 0.1878, 59.2903, 'Estay'}];

network_pairs{3} = [{0.3679, 0.2737, 30.0436/2, 'Eswitch'};
    {0.3161, 0.2482, 175, 'Estay'}];

network_pairs{4} = [{0.2800, 0.4228, 47.0955/2, 'Eswitch'};
    {0.2355, 0.4250, 175, 'Estay'}];

network_pairs{5} =  [{0.2921, 0.6927, 45/2, 'Eswitch'};
    {0.2119, 0.6799, 175, 'Estay'}];


%---------------------------------------
durations = cell(5,1);
%...
durations{1} = [{143.4682};
    {2.1808}];
durations{2} = [{187.2613};
    {2.0959}];

durations{3} = [{189.1735};
    {2.0199}];

durations{4} = [{165.4364};
    {2.1626}];

durations{5} =  [{179.6920};
    {2.2800}];

save(fullfile(sv_dir,'network_pairs.mat'),'network_pairs','durations')
