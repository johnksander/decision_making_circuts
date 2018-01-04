load('PS_parsweep_baseline_5.mat')
timestep = .25e-3; %.25 milisecond timestep
sim_results = sim_results{1}
sw = sim_results{2};
st = sim_results{1};
st_durs = cell2mat(st(:,1)) * timestep;
timestep = .25e-3; %.25 milisecond timestep
st_durs = cell2mat(st(:,1)) * timestep;
sw_durs = cell2mat(sw(:,1)) * timestep;
figure(2)
histogram(st_durs,80)
hold on
histogram(sw_durs,80)