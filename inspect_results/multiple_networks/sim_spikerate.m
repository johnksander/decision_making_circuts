function cell_raster = sim_spikerate(cell_raster,timestep,binsz)

%binsz = 2e-3;
num_binsamps = binsz./timestep; %num samples in Xms (the bin size)
if round(num_binsamps) - num_binsamps < 1e-12 %roundoff errors...
    num_binsamps = round(num_binsamps);
end
if mod(num_binsamps,1) ~= 0 %not evenly divisible... find next best thing 
    alt_binsz = binsz-(5e-3):1e-3:binsz+(5e-3);
    even_div = rem(alt_binsz, timestep) == 0;
    alt_binsz = alt_binsz(even_div);
    [~,binsz] = min(abs(binsz - alt_binsz)); %find closest evenly divisible binsize
    binsz = alt_binsz(binsz);
    num_binsamps = binsz/timestep;
end
 
raster_sz = size(cell_raster);
if mod(raster_sz(2),num_binsamps) ~= 0 %you have to trim it down, equally divisible by bin size
    cell_raster = cell_raster(:,1:end - mod(raster_sz(2),num_binsamps));
    raster_sz = size(cell_raster); %should be good now
end
bin_magic = [numel(cell_raster(:,1)), num_binsamps, numel(cell_raster(1,:))/num_binsamps]; %set up for a magic trick
cell_raster = reshape(cell_raster,bin_magic);
cell_raster = squeeze(sum(cell_raster,2)) ./ (num_binsamps * timestep); %convert to Hz
cell_raster = repmat(cell_raster,[num_binsamps 1 1]);
cell_raster = reshape(cell_raster,raster_sz); %put the rabit back in the hat
end
