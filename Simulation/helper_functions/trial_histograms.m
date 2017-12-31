function trial_histograms(options,trial_data,trialidx)

switch options.trial_hists
    case 'on'
        
        orient landscape
        fontsz = 16;
        for idx = 1:options.num_stims
            stim_current = options.trial_currents(trialidx,idx);
            ax(idx) =  subplot(options.num_stims,1,idx);
            histogram(trial_data{idx} * options.conv2secs)
            xlabel('state duration (s)','Fontsize',fontsz)
            ylabel('frequency','Fontsize',fontsz)
            tile4plot = {options.stim_labels{idx}, [' (' num2str(stim_current)...
               ' x base current)']};
            title(tile4plot,'Fontsize',fontsz)
            set(ax(idx),'fontsize',fontsz)
        end
        
        xlim = cell2mat(get(ax,'XLim'));
        xlim = [0 max(xlim(:))];
        ylim = cell2mat(get(ax,'YLim'));
        ylim = [0 max(ylim(:))];
        set(ax,'XLim',xlim)
        set(ax,'YLim',ylim)
        
        fig_name = ['trial_' num2str(trialidx)];
        fig_dir = fullfile(options.save_dir,[options.sim_name '_trialfigs']);
        if ~isdir(fig_dir)
            mkdir(fig_dir)
        end
        
        print(fullfile(fig_dir,fig_name),'-djpeg')
        pause(.1)
end