clear;clc;close all

%this simply removes the datapoints from spikerate figures 
%that aren't included in the the state duration parsweep data 

basedir = '~/Desktop/work/ACClab/rotation/project/';
res = 'rates_excit';
figdir = fullfile(basedir,'Results','figures_parsweep_D2t-slower_spikerates',res);


%this is the heatmap var from state durations parsweep
HM = load(fullfile(basedir,'Results','figures_parsweep_D2t_very_slow_baseline',...
    'logmean_duration','heatmap_nointerp.mat')); 
HM = HM.HM;

F = openfig(fullfile(figdir,'heatmap_nointerp.fig'));
FIM = get(gca,'Children');
FIM.AlphaData = ~isnan(HM);
FIM.CData(isnan(HM)) = NaN;

cb = colorbar;
cb.Label.String = 'Hz';

set(gcf, 'renderer', 'painters')
print(fullfile(figdir,'heatmap_nointerp_masked'),'-djpeg','-r400')
