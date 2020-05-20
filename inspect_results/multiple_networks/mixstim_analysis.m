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
figure; 
plt_order = [1,3,2,4];
matblue = [0,0.4470,0.7410];
matorange = [0.8500,0.3250,0.0980];
for idx = 1:num_types
    subplot(num_types,2,idx)
    curr_type = net_types{idx};
    if idx == 1
        col = matblue;
    elseif idx == 2
        col = matorange;
    end
    curr_data = data(ismember(data.type,curr_type),:);
    sc(idx) = scatter(curr_data.ratio,curr_data.duration,'MarkerEdgeColor',col);
    xlabel('ratio','FontWeight','b')
    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    set(gca,'FontSize',fz)
    subplot(num_types,2,idx+2) 
    scatter(curr_data.diff,curr_data.duration,'MarkerEdgeColor',col)
    xlabel('difference','FontWeight','b')
    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    set(gca,'FontSize',fz)
end
ax = findobj(gcf,'Type','axes');
linkaxes(ax(1:2),'y');linkaxes(ax(3:4),'y')
axP = get(gca,'Position');
[lp,l] = legend(sc,net_types,'FontWeight','b',...
    'Location','northoutside','Box','off','Orientation','horizontal','FontSize',fz+4);
set(gca, 'Position', axP)
lp.Position = [((1-lp.Position(3))/2),1-lp.Position(4),lp.Position(3:4)];

print(fullfile(datadir,'scatters'),'-djpeg')


figure; 
markers = {'x','+','^','v','d'};
sc = [];
legend_labs = {};
pidx = 0;
for idx = 1:num_types
    subplot(num_types,2,idx);hold on
    curr_type = net_types{idx};
    if idx == 1
        col = matblue;
    elseif idx == 2
        col = matorange;
    end
    curr_data = data(ismember(data.type,curr_type),:);
    sc(idx) = scatter(NaN,NaN,'MarkerEdgeColor',col,'Linewidth',1.25);
    
    Smags = unique(curr_data.total); %intensities
    for kidx = 1:numel(Smags)
        curr_mag = Smags(kidx);
        sm = curr_data.total == curr_mag;
        pidx = pidx + 1;
        scatter(curr_data.ratio(sm),curr_data.duration(sm),'MarkerEdgeColor',col,...
            'Marker',markers{kidx},'Linewidth',1.25);
    end

    xlabel('ratio','FontWeight','b')
    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    set(gca,'FontSize',fz)
    subplot(num_types,2,idx+2) ;hold on
    for kidx = 1:numel(Smags)
        sm = curr_data.total == Smags(kidx);
        pidx = pidx + 1;
        scatter(curr_data.diff(sm),curr_data.duration(sm),'MarkerEdgeColor',col,...
            'Marker',markers{kidx},'Linewidth',1.25);
    end
    xlabel('difference','FontWeight','b')
    if idx == 1,ylabel('log_{10}(s) sampling','FontWeight','b');end
    set(gca,'FontSize',fz)
end

ax = findobj(gcf,'Type','axes');
linkaxes(ax(1:2),'y');linkaxes(ax(3:4),'y')
axP = get(gca,'Position');
[lp,l] = legend(sc,net_types,'FontWeight','b',...
    'Location','northoutside','Box','off','Orientation','horizontal','FontSize',fz+4);
set(gca, 'Position', axP)
lp.Position = [((1-lp.Position(3))/2),1-lp.Position(4),lp.Position(3:4)];


print(fullfile(datadir,'scatters-2'),'-djpeg')



for idx = 1:num_types
    subplot(num_types,2,idx)
    curr_type = net_types{idx};
    fprintf('\n\n\n%s\n-----------\n',curr_type)
    curr_data = data(ismember(data.type,curr_type),:);
    curr_data(:,{'total','type'}) = []; %drop for stepwise lm
    stepwiselm(curr_data,'ResponseVar','duration','upper','quadratic')
    

end



%drop total from dataset for modeling...

curr_type ='fast';
curr_data = data(ismember(data.type,curr_type),:);
fitlm(curr_data,'duration ~ diff + ratio + diff^2 +  ratio^2')

curr_type ='slow';
curr_data = data(ismember(data.type,curr_type),:);
fitlm(curr_data,'duration ~ ratio + diff')

curr_type ='slow';
curr_data = data(ismember(data.type,curr_type),:);
fitlm(curr_data,'duration ~ diff^2')

curr_type ='slow';
curr_data = data(ismember(data.type,curr_type),:);
fitlm(curr_data,'duration ~ diff')


keyboard











