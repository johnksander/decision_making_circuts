function celltype = celltype_logicals(pool_options)
%make logicals for different simulation cell types
%inputs: number of total cells, proportion stay & switch out of total,
%proportion excitable & inhibitory out of each pool


num_stay = pool_options.num_cells * pool_options.sz_pools(1);
num_switch = pool_options.num_cells * pool_options.sz_pools(2);
num_excitable = pool_options.sz_EI(1) * num_stay;

cellvec = [1:pool_options.num_cells]';
celltype.pool_stay = cellvec <= num_stay; %stay pool
celltype.pool_switch = ~celltype.pool_stay; %switch pool

Estay = cellvec(1:num_stay) <= pool_options.sz_EI(1) * num_stay;
Eswitch = cellvec(num_stay+1:end) <= (pool_options.sz_EI(1) * num_switch) + num_stay;
celltype.excit = logical([Estay;Eswitch]); %excitatory
celltype.inhib = ~celltype.excit; %inhibitory

if abs(sum(pool_options.sz_EI + pool_options.sz_pools) - 2) > eps %this got weird with 2/3 & 1/3 proportions
    disp(sprintf('WARNING: pool size mismatch'))
end
