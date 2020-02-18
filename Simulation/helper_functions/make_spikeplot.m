function spikeplot = make_spikeplot(spikes)
%adds points to left/right of spike matrix, so spikes are easily visible
%with imagesc()

num_addspikes = 25; %how many spikes to add to left & right
num_addspikes = ones(1,num_addspikes);
kernel = [num_addspikes 1 num_addspikes]; %smoothing kernel
spikeplot = NaN(size(spikes));
for idx = 1:numel(spikes(:,1)) %for each cell
    spikeplot(idx,:) = conv(spikes(idx,:),kernel,'same');
end
spikeplot(spikeplot > 1) = 1; %just in case spikes were right next to each other or something 