function rateplot = make_spikerate_plot(spikes,cellpool,timestep)
%takes spike matrix % pool logical, 
%finds pool average 1/ISI for every timepoint

spikes = spikes(cellpool,:);
num_cells = numel(spikes(:,1));
cell_rates = NaN(size(spikes));

for idx = 1:num_cells
    spikecourse = spikes(idx,:);
    spikeinds = find(spikecourse);
    ISI = diff(spikeinds);
    ISI = 1 ./ (ISI .* timestep);
    for spikeidx = 1:numel(spikeinds)-1
        cell_rates(idx,spikeinds(spikeidx):spikeinds(spikeidx+1)) = ISI(spikeidx);
    end
end
rateplot = nanmean(cell_rates);












