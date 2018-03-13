clear
clc
format compact

load('sweep_params.mat')

sv_dir = 'helper_functions'; %save to a data structure..

data = [ItoE,EtoI,mu_duration];

fast_cuttoff = 1; %exclude networks switching under mean = 1 second
fast_cuttoff = mu_duration < fast_cuttoff;
fprintf('networks faster than min switch time = %i\n',sum(fast_cuttoff))
data = data(~fast_cuttoff,:);

data = sortrows(data,1); %sorted on ItoE


showme = data(:,2) < .72 & data(:,2) > .67;
sortrows(data(showme,:),3)



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
