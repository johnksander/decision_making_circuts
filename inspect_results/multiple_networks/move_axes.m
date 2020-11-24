clear; clc
format compact; close all

%this takes one network pair out of the larger figure and replots by itself

basedir = '~/Desktop/work/ACClab/rotation/project';
redir = fullfile(basedir,'Results/figures_nets_mixstim/durations/Nmin_250');
fns = {'nets_mixstim_total_time_log.fig','nets_mixstim_total_time_log_diff.fig'};

fz = 20; %12;
Ntick = 5;

fast_inds = 13:14;
slow_inds = 15:16;
ptime = 1;

f2 = figure(2); set(gcf,'Renderer','painters');
scr = get( groot,'Screensize');
fwid = 1000;
pos = [1,scr(3),fwid,fwid*.4];
set(gcf,'Position',pos);

Np = 2; sax = [];
for idx = 1:Np
    sax(idx) = subplot(1,2,idx);
end


%ratio, difference
for fidx = 1:numel(fns)
    
    f1 = openfig(fullfile(redir,fns{fidx}));
    
    ch1 = get(f1,'Children');
    %only get fast network from each figure 
    idx = fidx + 0;  %previously: idx = ((fidx-1)*numel(fns)) + 1;
    
    figure(f2);set(gcf,'Renderer','painters')
    oldax = ch1(fast_inds);
    oldax(2).Position = get(sax(idx),'Position');pause(ptime)
    leglab = oldax(1).String;
    oldax = oldax(2);
    newax = copyobj(oldax,f2);
    delete(sax(idx))
    figure(f2);axes(newax);
    
    %newax.PositionConstraint = 'outerposition';
    hold on
         
    legend(leglab,'Location','northeast','Box','off');
    %if idx == 1 %No need to label these anymore
    %    title('fast network','FontWeight','b','FontSize',fz+2)
    %end
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',fz)
    ylabel('seconds (log scale)','FontWeight','b')
    if fidx == 1
        Xtick = get(gca,'XTickLabel');
        curr_Xtick = get(gca,'XTick');
        if numel(curr_Xtick) ~= numel(oldax.XTick) %errgghh
            set(gca,'XTick',oldax.XTick)
            set(gca,'XTickLabel',Xtick)
        end
        xlabel('aversive / hedonic','FontWeight','b')
    else
        xlabel('aversive - hedonic','FontWeight','b')
    end
    
    Ytick = get(gca,'YTick');
    Ytick = linspace(Ytick(1),Ytick(end),Ntick);
    set(gca,'YTick',Ytick)
    Ytick = 10.^Ytick;
    Ytick = cellfun(@(x) sprintf('%.1f',x),num2cell(Ytick),'UniformOutput',false);
    Ytick = regexprep(Ytick,'.0$','');
    Ytick = regexprep(Ytick,'^0.','.');
    set(gca,'YTickLabel',Ytick);
    
    
    pos = get(gca,'Position');
    %set(gca,'Position',[pos(1)*.8,pos(2)*1.2,pos(3:4)]);
    set(gca,'Position',[pos(1),pos(2)*1.5,pos(3:4)]);
    
    
    
    close(f1)
end

print('fig7-stim_mixture','-djpeg','-r400')



%previously, we did a four-paneled plot showing the slow and fast networks
% 
% %ratio, difference
% for fidx = 1:numel(fns)
%     
%     f1 = openfig(fullfile(redir,fns{fidx}));
%     
%     ch1 = get(f1,'Children');
%     %fast left, slow right
%     idx = ((fidx-1)*numel(fns)) + 1;
%     
%     figure(f2);set(gcf,'Renderer','painters')
%     oldax = ch1(fast_inds);
%     oldax(2).Position = get(sax(idx),'Position');pause(ptime)
%     leglab = oldax(1).String;
%     oldax = oldax(2);
%     newax = copyobj(oldax,f2);
%     delete(sax(idx))
%     figure(f2);axes(newax)
%     set(gca,'FontSize',fz)
%     legend(leglab,'Location','northeast','Box','off','FontSize',fz-2);
%     Ytick = get(gca,'YTick');
%     Ytick = linspace(Ytick(1),Ytick(end),Ntick);
%     set(gca,'YTick',Ytick)
%     Ytick = 10.^Ytick;
%     Ytick = cellfun(@(x) sprintf('%.1f',x),num2cell(Ytick),'UniformOutput',false);
%     Ytick = regexprep(Ytick,'.0$','');
%     Ytick = regexprep(Ytick,'^0.','.');
%     set(gca,'YTickLabel',Ytick);
%     ylabel('seconds (log scale)')
%     if idx == 1
%         title('fast network','FontWeight','b','FontSize',fz+2)
%     end
%     set(gca,'TickLabelInterpreter','latex')
%     if fidx == 1
%         Xtick = get(gca,'XTickLabel');
%         curr_Xtick = get(gca,'XTick');
%         if numel(curr_Xtick) ~= numel(oldax.XTick) %errgghh
%             set(gca,'XTick',oldax.XTick)
%             set(gca,'XTickLabel',Xtick)
%         end
%         xlabel('aversive / hedonic')
%     else
%         xlabel('aversive - hedonic')
%     end
%     %text(-.1,1.1,char(idx + 64),'Units','normalized',...
%     %    'FontWeight','b','FontSize',fz)
%     
%     
%     
%     
%     idx = idx + 1;
%     oldax = ch1(slow_inds);
%     oldax(2).Position = get(sax(idx),'Position');pause(ptime)
%     leglab = oldax(1).String;
%     oldax = oldax(2);
%     newax = copyobj(oldax,f2);
%     delete(sax(idx))
%     figure(f2);axes(newax)
%     set(gca,'FontSize',fz)
%     legend(leglab,'Location','southwest','Box','off','FontSize',fz-2);
%     Ytick = get(gca,'YTick');
%     Ytick = linspace(Ytick(1),Ytick(end),Ntick);
%     set(gca,'YTick',Ytick)
%     Ytick = 10.^Ytick;
%     Ytick = cellfun(@(x) sprintf('%.1f',x),num2cell(Ytick),'UniformOutput',false);
%     Ytick = regexprep(Ytick,'.0$','');
%     Ytick = regexprep(Ytick,'^0.','.');
%     set(gca,'YTickLabel',Ytick);
%     ylabel('seconds (log scale)')
%     if idx == 2
%         title('slow network','FontWeight','b','FontSize',fz+2)
%     end
%     
%     set(gca,'TickLabelInterpreter','latex')
%     %setting the interpreter to latex seems to default "x/y" tick labels to
%     %the nicefrac{} latex format, so that's nice.
%     if fidx == 1
%         Xtick = get(gca,'XTickLabel');
%         curr_Xtick = get(gca,'XTick');
%         if numel(curr_Xtick) ~= numel(oldax.XTick) %errgghh
%             set(gca,'XTick',oldax.XTick)
%             set(gca,'XTickLabel',Xtick)
%         end
%         xlabel('aversive / hedonic')
%     else
%         xlabel('aversive - hedonic')
%     end
%     
%     %text(-.1,1.1,char(idx + 64),'Units','normalized',...
%     %    'FontWeight','b','FontSize',fz)
%     
%     close(f1)
% end
% 
% %set(gcf,'Renderer','painters')
% print('fig7-stim_mixture','-djpeg','-r400')





%this is the idea...
%https://www.mathworks.com/matlabcentral/answers/101273-how-can-i-put-existing-figures-in-different-subplots-in-another-figure-in-matlab-6-5-r13



% Xtick = get(gca,'XTickLabel');
% Xtick = strrep(Xtick,'.','');
% Xtick = strrep(Xtick,'/',':');
% set(gca, 'XTicklabel',Xtick);
% xlabel('aversive : hedonic')
% curr_Xtick = get(gca,'XTick');
% if numel(curr_Xtick) ~= numel(oldax.XTick) %errgghh
%     set(gca,'XTick',oldax.XTick)
%     set(gca,'XTickLabel',Xtick)
% end



%set(gca,'TickLabelInterpreter','latex')
%Xtick = get(gca,'XTickLabel');
%Xtick = strrep(Xtick,'.','');
%Xtick = cellfun(@(x) strsplit(x,'/'),Xtick,'UniformOutput',false);
%Xtick = cellfun(@(x) sprintf('$\\nicefrac{%s}{%s}$',x{1},x{2}),Xtick,'UniformOutput',false);




