function avg_rates = avg_poolrate(avg_rates,last_spike,spiking_cells,celltype,idx)
%get inds of last spikes, output average pool rate 

avg_rates.ISIs(spiking_cells) = 1 ./ last_spike(spiking_cells);


%this is kinda dumb but w/e, doesn't look like it has errors..
Estay = celltype.excit & celltype.pool_stay;
Eswitch = celltype.excit & celltype.pool_switch;
Istay = celltype.inhib & celltype.pool_stay;
Iswitch = celltype.inhib & celltype.pool_switch;

avg_rates.Estay(idx) = nanmean(avg_rates.ISIs(Estay));
avg_rates.Eswitch(idx) = nanmean(avg_rates.ISIs(Eswitch));
avg_rates.Istay(idx) = nanmean(avg_rates.ISIs(Istay));
avg_rates.Iswitch(idx) = nanmean(avg_rates.ISIs(Iswitch));
