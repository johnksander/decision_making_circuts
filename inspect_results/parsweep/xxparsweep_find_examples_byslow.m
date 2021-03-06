clear
clc
format compact
close all

%this matches exampes by slow networks

%specify simulation
%---sim setup-----------------
sim_name = 'parsweep_slowD_Rlim_baseline';
Flabel = 'slowD'; %network pairs will be saved with this filename 
basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project/';
figdir = fullfile(basedir,'Results',['figures_' sim_name]);
resdir = fullfile(basedir,'Results',sim_name);
helper_dir = fullfile(basedir,'helper_functions');
addpath(helper_dir)
svdir = fullfile(helper_dir,'network_pairs');
if ~isdir(svdir),mkdir(svdir);end
FN = fullfile(svdir,Flabel);

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
slow_range = log10([55,65]); %seconds
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

scatter(netX(cand_nets),netY(cand_nets),'black')
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



%test inds
param_perc = @(x,y) (range(y)*x )+ min(y); %return the Xth percentile value from parameter range Y


Npairs = 5;
network_pairs = cell(Npairs,1);

%find the one at the "top of the curve"
minmax_norm = @(x) (x - min(x)) ./ (max(x) - min(x));

center_pair = array2table(NaN(2,3),'VariableNames',{'ItoE','EtoI','duration'},'RowNames',{'slow','fast'});

XYcoords = [ItoE,EtoI];
XYcoords = num2cell(XYcoords,1);
XYcoords = cellfun(minmax_norm ,XYcoords,'UniformOutput',false);

slowP = [ItoE(slow_nets),EtoI(slow_nets)];
slowXY = cellfun(@(x) x(slow_nets),XYcoords,'UniformOutput',false);
slowXY = cat(2,slowXY{:});

fastP = [ItoE(fast_nets),EtoI(fast_nets)];
fastXY = cellfun(@(x) x(fast_nets),XYcoords,'UniformOutput',false);
fastXY = cat(2,fastXY{:});

minP = pdist2([0,0],slowXY)';
minP = find(minP == min(minP)); %closest to origin 
slow_coords = slowXY(minP,:); %save this for next step 
slowP = slowP(minP,:);
curr_net = ItoE == slowP(1) & EtoI == slowP(2);
curr_net = find(curr_net);
center_pair{'slow','ItoE'} = ItoE(curr_net);
center_pair{'slow','EtoI'} = EtoI(curr_net);
center_pair{'slow','duration'} = Sdurs(curr_net);

%minP = pdist2([0,0],fastXY)'; %can also do origin like above 
minP = pdist2(slow_coords,fastXY)'; %can also do origin like above 
minP = find(minP == min(minP)); %closest to slow net
fastP = fastP(minP,:);
curr_net = ItoE == fastP(1) & EtoI == fastP(2);
curr_net = find(curr_net);
center_pair{'fast','ItoE'} = ItoE(curr_net);
center_pair{'fast','EtoI'} = EtoI(curr_net);
center_pair{'fast','duration'} = Sdurs(curr_net);

network_pairs{5} = center_pair;


%find the rest
p.ItoE = ItoE;
p.EtoI = EtoI;
range_vals = [1, .5]; %get one pair at the extreme, one at half-way
p_types = {'ItoE','EtoI'};

pair_ind = 0;
for idx = 1:numel(p_types)
    
    curr_Ps = p.(p_types{idx}); %connection strength values for ItoE or EtoI
    switch p_types{idx}
        case 'ItoE'
            this_arm.fast = p.ItoE > center_pair{'fast','ItoE'};
            this_arm.slow = p.ItoE > center_pair{'slow','ItoE'} & p.EtoI < center_pair{'slow','EtoI'};
            %this_arm = p.ItoE > center_pair{'slow','ItoE'} & p.EtoI < center_pair{'slow','EtoI'};
            
        case 'EtoI'
            this_arm.fast = p.EtoI > center_pair{'fast','EtoI'};
            this_arm.slow = p.EtoI > center_pair{'slow','EtoI'} & p.EtoI > center_pair{'slow','EtoI'};
            %this_arm = p.ItoE < center_pair{'slow','ItoE'} & p.EtoI > center_pair{'slow','EtoI'};
            
    end
    
    for RVidx = 1:numel(range_vals)
        
        pair_ind = pair_ind + 1;
        net_pair = array2table(NaN(2,3),'VariableNames',{'ItoE','EtoI','duration'},'RowNames',{'slow','fast'});
        
        %find the slow network 
        curr_RV = range_vals(RVidx); %where in the parameter space 
        param_range = curr_Ps(slow_nets & this_arm.slow); %find the range for this "arm"
        
        curr_net = param_perc(curr_RV,param_range);
        curr_net = nearest_val(curr_net,param_range);
        curr_net = find(curr_Ps == curr_net & slow_nets);

        net_pair{'slow','ItoE'} = ItoE(curr_net);
        net_pair{'slow','EtoI'} = EtoI(curr_net);
        net_pair{'slow','duration'} = Sdurs(curr_net);
        
        %find the matching fast one
        curr_net = nearest_val(net_pair{'slow',p_types{idx}},curr_Ps(fast_nets & this_arm.fast));
        curr_net = find(curr_Ps == curr_net & fast_nets);
        
        net_pair{'fast','ItoE'} = ItoE(curr_net);
        net_pair{'fast','EtoI'} = EtoI(curr_net);
        net_pair{'fast','duration'} = Sdurs(curr_net);
        
        network_pairs{pair_ind} = net_pair;
    end
end

%mark pairs on the big figure 

X = cellfun(@(x) x{:,'ItoE'},network_pairs,'UniformOutput',false);
Y = cellfun(@(x) x{:,'EtoI'},network_pairs,'UniformOutput',false);
Z = cellfun(@(x) log10(x{:,'duration'}),network_pairs,'UniformOutput',false);
X = cat(1,X{:});Y = cat(1,Y{:});Z = cat(1,Z{:});

figure(SP) %return focus to this
hold on

scatter3(X,Y,Z,1000,'red','p','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1); %mark 'em w/ stars
hold off

  
%now for the heatmap

%need a logical showing all the networks you picked... 
HMcoords = [ItoE,EtoI];
HMcoords = ismember(HMcoords,[X,Y],'rows');

figure(HM)
hold on

scatter(netX(HMcoords),netY(HMcoords),500,'red','p','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
hold off


%save the network pair data
save(FN,'network_pairs')


% net_pair = array2table(NaN(2,3),'VariableNames',{'ItoE','EtoI','duration'},'RowNames',{'slow','fast'});
%
% curr_net = param_perc(1,ItoE(slow_nets));
% curr_net = nearest_ind(curr_net,ItoE);
%
% net_pair{'slow','ItoE'} = ItoE(curr_net);
% net_pair{'slow','EtoI'} = EtoI(curr_net);
% net_pair{'slow','duration'} = Sdurs(curr_net);
%
% %find the matching fast one
% curr_net = nearest_val(net_pair{'slow','ItoE'},ItoE(fast_nets));
% curr_net = find(ItoE == curr_net & fast_nets);
%
% net_pair{'fast','ItoE'} = ItoE(curr_net);
% net_pair{'fast','EtoI'} = EtoI(curr_net);
% net_pair{'fast','duration'} = Sdurs(curr_net);


% scatter3(ItoE(cand_nets),EtoI(cand_nets),log10(state_dur(cand_nets)),20,'black','filled',...
%     'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);pause(1) %mark 'em


%plot3(ItoE,EtoI,outcome,'.','MarkerSize',15,'color','blue') %give them outlines
%plot3(ItoE,EtoI,outcome,'.','MarkerSize',10,'color','green')

%
% scatter3(Nplots(:,1),Nplots(:,2),Nplots(:,3),1000,'red','p','filled',...
%     'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1); %mark 'em w/ stars



%from the last time around::

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
% 
% %original --------------------
% network_pairs = cell(5,1);
% %...
% network_pairs{1} = [{1.2904, 0.1948, 75, 'Eswitch'};
%     {1.2877,  0.1734, 100, 'Estay'}];
% 
% network_pairs{2} = [{0.8069, 0.2067, 49.7632, 'Eswitch'};
%     {0.7936, 0.1878, 59.2903, 'Estay'}];
% 
% network_pairs{3} = [{0.3679, 0.2737, 30.0436, 'Eswitch'};
%     {0.3161, 0.2482, 175, 'Estay'}];
% 
% network_pairs{4} = [{0.2800, 0.4228, 47.0955, 'Eswitch'};
%     {0.2355, 0.4250, 175, 'Estay'}];
% 
% network_pairs{5} =  [{0.2921, 0.6927, 45, 'Eswitch'};
%     {0.2119, 0.6799, 175, 'Estay'}];
% %end originial --------------------
% 
% %now trying to equate------------------
% %03212018: I'm just halfing each E-switch stim
% network_pairs{1} = [{1.2904, 0.1948, 75/2, 'Eswitch'};
%     {1.2877,  0.1734, 100, 'Estay'}];
% 
% network_pairs{2} = [{0.8069, 0.2067, 49.7632/2, 'Eswitch'};
%     {0.7936, 0.1878, 59.2903, 'Estay'}];
% 
% network_pairs{3} = [{0.3679, 0.2737, 30.0436/2, 'Eswitch'};
%     {0.3161, 0.2482, 175, 'Estay'}];
% 
% network_pairs{4} = [{0.2800, 0.4228, 47.0955/2, 'Eswitch'};
%     {0.2355, 0.4250, 175, 'Estay'}];
% 
% network_pairs{5} =  [{0.2921, 0.6927, 45/2, 'Eswitch'};
%     {0.2119, 0.6799, 175, 'Estay'}];
% 
% 
% %---------------------------------------
% durations = cell(5,1);
% %...
% durations{1} = [{143.4682};
%     {2.1808}];
% durations{2} = [{187.2613};
%     {2.0959}];
% 
% durations{3} = [{189.1735};
%     {2.0199}];
% 
% durations{4} = [{165.4364};
%     {2.1626}];
% 
% durations{5} =  [{179.6920};
%     {2.2800}];
% 
% save(fullfile(sv_dir,'network_pairs.mat'),'network_pairs','durations')
