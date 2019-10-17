function [options,outcome] = check_rate_limit(spikes,celltype,options)
%this performs ratelimit check, and checks for oscillators that passed
%initial bistability check

outcome.status = 'pass';
if ~strcmp(options.record_spiking,'ratelim_only')
    %if we're recording spikes for the whole experiment, must select ratelim window.
    %"ratelim_only" now defaults to only recording during the window.
    check_start = options.ratelim.start / options.timestep;
    check_stop = options.ratelim.stop / options.timestep;
    spikes = spikes(:,check_start:check_stop);
end

%get rates
window_sz = 75e-3;
raster_sz = size(spikes);
num_binsamps = window_sz/options.timestep; %num samples in 2ms
if round(num_binsamps) - num_binsamps < 1e-12 %roundoff errors...
    num_binsamps = round(num_binsamps);
end
if mod(raster_sz(2),num_binsamps) ~= 0 %you have to trim it down, equally divisible by bin size
    Ntrim = mod(raster_sz(2),num_binsamps);
    spikes = spikes(:,1+Ntrim:end);
    raster_sz = size(spikes); %should be good now
end
bin_magic = [raster_sz(1), num_binsamps,raster_sz(2)/num_binsamps]; %set up for a magic trick
spikes = reshape(spikes,bin_magic);
spikes = squeeze(sum(spikes,2)) ./ (num_binsamps * options.timestep); %convert to Hz
spikes = repmat(spikes,[num_binsamps 1 1]);
spikes = reshape(spikes,raster_sz); %put the rabit back in the hat

%normal people indexing that makes sense, then take the mean
rate.Estay = mean(spikes(celltype.excit & celltype.pool_stay,:),1);
rate.Eswitch = mean(spikes(celltype.excit & celltype.pool_switch,:),1);
rate.Istay = mean(spikes(celltype.inhib & celltype.pool_stay,:),1);
rate.Iswitch = mean(spikes(celltype.inhib & celltype.pool_switch,:),1);

Imax = max([rate.Istay;rate.Iswitch]);
Emax = max([rate.Estay;rate.Eswitch]);

Eover = Emax > options.ratelim.E;
Eover = overlimit_times(Eover,options.timestep);

Iover = Imax > options.ratelim.I;
Iover = overlimit_times(Iover,options.timestep);

%record approximate rates
outcome.Erate = mean(Emax);
outcome.Irate = mean(Imax);


if sum(Eover > options.ratelim.tmax) > 0
    update_logfile(':::Sustained E-spiking over limit:::',options.output_log)
    message = sprintf('---limit: %iHz for %.2fs',options.ratelim.E,options.ratelim.tmax);
    update_logfile(message,options.output_log)
    outcome.status = 'fail';
    
elseif outcome.Erate > options.ratelim.E * options.ratelim.mulim
    update_logfile(':::mean E-spiking over limit:::',options.output_log)
    message = sprintf('---mean limit: %.2fHz',options.ratelim.E * options.ratelim.mulim);
    update_logfile(message,options.output_log)
    update_logfile(sprintf('---mean E-rate: %.2fHz',outcome.Erate),options.output_log)
    outcome.status = 'fail';
end


if sum(Iover > options.ratelim.tmax) > 0
    update_logfile(':::Sustained I-spiking over limit:::',options.output_log)
    message = sprintf('---limit: %iHz for %.2fs',options.ratelim.I,options.ratelim.tmax);
    update_logfile(message,options.output_log)
    outcome.status = 'fail';
    
elseif outcome.Irate > options.ratelim.I * options.ratelim.mulim
    update_logfile(':::mean I-spiking over limit:::',options.output_log)
    message = sprintf('---mean limit: %.2fHz',options.ratelim.I * options.ratelim.mulim);
    update_logfile(message,options.output_log)
    update_logfile(sprintf('---mean I-rate: %.2fHz',outcome.Irate),options.output_log)
    outcome.status = 'fail';
end


%reset options.record_spiking to off, if ratelim_only
if strcmp(options.record_spiking,'ratelim_only')
    options.record_spiking = 'off';
end

    function v = overlimit_times(x,t)
        %takes---
        %logical vector x : rates over limit
        %t: timestep
        %returns vector of consequtive timepoints over limit
        q = diff([0 x 0]);
        v = find(q == -1) - find(q == 1);
        v = v * t;
    end
end


%for plotting & checking this functionality out

% figure;
% lnsz = 3;fontsz = 12;
% plot(rate.Estay,'Linewidth',lnsz);hold on
% plot(rate.Eswitch,'Linewidth',lnsz)
% plot(rate.Istay,'Linewidth',lnsz)
% plot(rate.Iswitch,'Linewidth',lnsz)
% Xticks = num2cell(get(gca,'Xtick'));title('Spiking')
% Xlabs = cellfun(@(x) sprintf('%.1f',x*options.timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)});xlabel('time (s)')
% legend({'E-stay','E-switch','I-stay','I-switch'},'location','northoutside','Orientation','horizontal')
% set(gca,'FontSize',fontsz);axis tight;hold off
% figure;
% plot(Emax,'Linewidth',lnsz);hold on
% plot(Imax,'Linewidth',lnsz)
% Xticks = num2cell(get(gca,'Xtick'));title('Spiking')
% Xlabs = cellfun(@(x) sprintf('%.1f',x*options.timestep),Xticks,'UniformOutput', false); %this is for normal stuff
% set(gca,'Xdir','normal','Xtick',cell2mat(Xticks),'XTickLabel', Xlabs);
% ylabel({sprintf('Mean pool Hz  (%ims bins)',window_sz*1e3)});xlabel('time (s)')
% legend({'E-max','I-max'},'location','northoutside','Orientation','horizontal')
% set(gca,'FontSize',fontsz);axis tight;hold off

% %this wont give static workplace error
% plot(rate.Estay,'Linewidth',3);hold on
% plot(rate.Eswitch,'Linewidth',3)
% plot(rate.Istay,'Linewidth',3)
% plot(rate.Iswitch,'Linewidth',3)
% set(gca,'Xdir','normal','Xtick',cell2mat(num2cell(get(gca,'Xtick'))),...
%     'XTickLabel',cellfun(@(x) sprintf('%.1f',x*options.timestep),...
%     num2cell(get(gca,'Xtick')),'UniformOutput', false))
% figure;
% plot(Emax,'Linewidth',3);hold on
% plot(Imax,'Linewidth',3)
% set(gca,'Xdir','normal','Xtick',cell2mat(num2cell(get(gca,'Xtick'))),...
%     'XTickLabel',cellfun(@(x) sprintf('%.1f',x*options.timestep),...
%     num2cell(get(gca,'Xtick')),'UniformOutput', false))
