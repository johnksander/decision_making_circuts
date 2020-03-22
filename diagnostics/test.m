
clear
clc
format compact
hold off;close all
%investigating model behavior


addpath('../')
options = struct();
options.netpair_file = 'D2t-slower';

stims = array2table(NaN(5,2),'VariableNames',{'Estay','Eswitch'});

for idx = 1:10
    
    net = get_network_params(idx,options);
    
    stims.(net.stim_targs{1})(ceil(idx./2)) = net.trial_stimuli{1}(1);
end
keyboard
%do 3 & 4


scatter(stims.Estay,stims.Eswitch,100,'filled')