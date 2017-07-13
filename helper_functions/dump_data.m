function dump_data(rolling_timecourse,celltype,rtSpikes,rtNoise,dumpfn)
%dump switch timecourse data to a text file
%spike & noise data
%judges' score: -50 style points for this function. 


dumpdata = NaN(8,numel(rolling_timecourse(1,:,1))); %one row for each kinda cell (x2: spikes, then noise)
Estay = celltype.excit & celltype.pool_stay;
Istay = celltype.inhib & celltype.pool_stay;
Eswitch = celltype.excit & celltype.pool_switch;
Iswitch = celltype.inhib & celltype.pool_switch;

dumpdata(1,:) = sum(rolling_timecourse(Estay,:,rtSpikes)); %this is stupid but just go with it
dumpdata(2,:) = sum(rolling_timecourse(Istay,:,rtSpikes));
dumpdata(3,:) = sum(rolling_timecourse(Eswitch,:,rtSpikes));
dumpdata(4,:) = sum(rolling_timecourse(Iswitch,:,rtSpikes));

dumpdata(5,:) = sum(rolling_timecourse(Estay,:,rtNoise)); %again... just go with it.
dumpdata(6,:) = sum(rolling_timecourse(Istay,:,rtNoise));
dumpdata(7,:) = sum(rolling_timecourse(Eswitch,:,rtNoise));
dumpdata(8,:) = sum(rolling_timecourse(Iswitch,:,rtNoise));

dlmwrite(dumpfn,dumpdata);

