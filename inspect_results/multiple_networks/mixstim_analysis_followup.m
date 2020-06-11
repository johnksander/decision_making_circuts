clear
clc
format compact
hold off;close all


basedir = '~/Desktop/work/ACClab/rotation/project'; %'~/Desktop/work/ACClab/rotation/project/';
datadir = fullfile(basedir,['Results/' ...
    'figures_nets_mixstim/durations/Nmin_250/analysis-logmu']);
%'figures_nets_mixstim_netpair-2/durations/Nmin_250/analysis-logmu']);

fz = 16;

data = load(fullfile(datadir,'model_data'));
data = data.Mdata;
net_types = {'slow','fast'}; %unique(data.type);
num_types = numel(net_types);
mixes = {'diff','ratio'};
plt_order = [1,3,2,4];
matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
matyellow = [0.9290,.6940,.1250];

data = data(data.net_index == 2,:);

% data = data(ismember(data.type,'fast'),:);
% data = data(ismember(round(data.total),[198,396,793]),:);

figure;
orient landscape %single pair 
lnsz = 1;
lntype = '-.';
sc = [];
legend_labs = {};
pidx = 0;
for idx = 1:numel(mixes)
    subplot(1,numel(mixes),idx);hold on
    
    curr_mix = mixes{idx};
    curr_type = net_types{idx};
    if idx == 1
        col = matblue;
    elseif idx == 2
        col = matorange;
    end
    %curr_data = data(ismember(data.type,curr_type),:);
    
    curr_data = data;
    Smags = unique(curr_data.total); %intensities
    for kidx = 1:numel(Smags)
        curr_mag = Smags(kidx);
        sm = curr_data.total == curr_mag;
        switch curr_mix
            case 'ratio'
                plot(curr_data.ratio(sm),curr_data.duration(sm),'LineWidth',3)
            case 'diff'
                plot(curr_data.diff(sm),curr_data.duration(sm),'LineWidth',3)
        end
    end
    
    ylim([-1,2])

    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    xlabel(curr_mix,'FontWeight','b')
    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    set(gca,'FontSize',fz)
end


close all

Smags = unique(curr_data.total); %intensities
cols = [matblue;matorange;matyellow];
figure;scatter3(NaN,NaN,NaN);hold on
pl = [];
for kidx = 1:numel(Smags)
    curr_mag = Smags(kidx);
    sm = curr_data.total == curr_mag;
    z = curr_data.duration(sm);
    y = curr_data.ratio(sm);
    x = curr_data.diff(sm);
    c = cols(kidx,:);
    pl(kidx) = scatter3(x,y,z,'filled','MarkerEdgeColor',c,'MarkerFaceColor',c);
end

zlabel('log_{10}(s) sampling','FontWeight','bold')
ylabel('E-switch / E-stay','FontWeight','bold')
xlabel( 'E-switch - E-stay (Hz)','FontWeight','bold')

ylim([0,1])
yticks('manual')
ytick = num2cell(get(gca,'YTick'));
ytick = cellfun(@(x) sprintf('%1.1f/%.1f',x,1-x),ytick,'UniformOutput',false);
ytick = strrep(ytick,'0.','.');
ytick = strrep(ytick,'1.0','1');
ytick = strrep(ytick,'.0','0');
set(gca,'YTickLabel',ytick)


legend_labs = cellfun(@(x) sprintf('%.0f Hz total',x),num2cell(Smags),'UniformOutput',false);
pause(1);legend(pl,legend_labs,'Location','best','Box','off');pause(1)


title(sprintf('network #%i\n%s',curr_data.pair_index),'FontWeight','bold')
keyboard
view(90,0)
view(360,0)
        
base_targ_cells = 'Eswitch'; %base everything off this
stimtarg_vals = {'baseline','Estay','Eswitch'}; %this is dumb
stimtarg_labels = {'baseline','fast','slow'};
base_targ_labels = stimtarg_labels(strcmp(stimtarg_vals,base_targ_cells));
      
if plt_idx == num_net_types-1 || plt_idx == num_net_types
    mix_info = {'Estay','Eswitch'};
    mix_info = sprintf('%s / %s',mix_info{strcmp(mix_info,base_targ_cells)},...
        mix_info{~strcmp(mix_info,base_targ_cells)});
    mix_info = strrep(mix_info,'E','E-');
    switch opt.Xax
        case 'diff'
            mix_info = strrep(mix_info,' / ',' - ');
            mix_info = [mix_info ' (Hz)'];
    end
    xlabel(mix_info,'FontWeight','bold')
    
end

        
        ytick = num2cell(get(gca,'XTick'));
        ytick = cellfun(@(x) sprintf('%1.1f/%.1f',x,1-x),ytick,'UniformOutput',false);
        ytick = strrep(ytick,'0.','.');
        ytick = strrep(ytick,'1.0','1');
        ytick = strrep(ytick,'.0','0');
        set(gca,'XTickLabel',ytick)
        
switch opt.Xax
    case 'ratio'
        axis tight;xlim([0:1])
        legend_labs = cellfun(@(x) sprintf('%.0f Hz total',x),num2cell(Smags),'UniformOutput',false);
        pause(1);legend(legend_labs,'Location','best','Box','off');pause(1)
        ytick = num2cell(get(gca,'XTick'));
        ytick = cellfun(@(x) sprintf('%1.1f/%.1f',x,1-x),ytick,'UniformOutput',false);
        ytick = strrep(ytick,'0.','.');
        ytick = strrep(ytick,'1.0','1');
        ytick = strrep(ytick,'.0','0');
        set(gca,'XTickLabel',ytick)
    case 'diff'
        axis tight;
        legend_labs = cellfun(@(x) sprintf('%.0f Hz total',x),num2cell(Smags),'UniformOutput',false);
        pause(1);legend(legend_labs,'Location','best','Box','off');pause(1)
end

subplot(1,numel(mixes),idx);hold on

curr_mix = mixes{idx};
curr_type = net_types{idx};
if idx == 1
    col = matblue;
elseif idx == 2
    col = matorange;
end
%curr_data = data(ismember(data.type,curr_type),:);

curr_data = data;
Smags = unique(curr_data.total); %intensities
for kidx = 1:numel(Smags)
    curr_mag = Smags(kidx);
    sm = curr_data.total == curr_mag;
    switch curr_mix
        case 'ratio'
            plot(curr_data.ratio(sm),curr_data.duration(sm),'LineWidth',3)
        case 'diff'
            plot(curr_data.diff(sm),curr_data.duration(sm),'LineWidth',3)
    end
end


if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
xlabel(curr_mix,'FontWeight','b')
if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
set(gca,'FontSize',fz)




keyboard

return









