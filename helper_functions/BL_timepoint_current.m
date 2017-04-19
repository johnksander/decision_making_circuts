function Iapp = BL_timepoint_current(options,current_info,state,durations,celltype)
%intended as a temporary function for baseline test....
%inputs: options, current, state, celltype parameter structures & timepoint index

%initpulse_time = 200; %idx is just index number, should be invariant across units (i.e. PMvsJK)
Iapp = current_info.noise_sigma * randn(current_info.num_cells,1);
Estay = celltype.pool_stay & celltype.excit;
Eswitch = celltype.pool_switch & celltype.excit;

switch options.modeltype
    case 'PM' %not intended for PM model
        
        %         %inhibitory neruons always get the same currrent here!!
        %         Iapp(celltype.inhib) = Iapp(celltype.inhib) + current_info.basecurr(celltype.inhib);
        %
        %         if current_info.timeidx <= current_info.initpulse_time %initial pulse window, force first switch
        %             %stay excitable cells get normal base current
        %             Iapp(Estay) = Iapp(Estay) + current_info.basecurr(Estay);
        %             %switch exctiables get half current, force first switch
        %             Iapp(Eswitch) = Iapp(Eswitch) + (.5 * current_info.basecurr(Eswitch));
        %
        %         else %past initial pulse window
        %             if sum(state.switch == state.now) == numel(state.now) %we're in a switch state
        %                 %stay excitables get base current, no stim bias
        %                 Iapp(Estay) = Iapp(Estay) + current_info.basecurr(Estay);
        %                 %switch excitables get no current
        %                 %Iapp(Eswitch) = Iapp(Eswitch) + 0;%...
        %             elseif sum(state.stay == state.now) == numel(state.now) %we're in a stay state
        %                 %switch excitables get full current, no matter which stim
        %                 Iapp(Eswitch) = Iapp(Eswitch) + current_info.basecurr(Eswitch);
        %                 %add stimulus-specific current
        %                 if mod(numel(durations{state.stay}),2) == 0 %stim A (count is one behind, based off already recorded duration..)
        %                     %stay excitables get full current adjusted by stim modifier
        %                     Iapp(Estay) = Iapp(Estay) + (current_info.basecurr(Estay) * current_info.stimA);
        %                 else %stim B
        %                     Iapp(Estay) = Iapp(Estay) + (current_info.basecurr(Estay) * current_info.stimB);
        %                 end
        %             end
        %         end
        %
        %         %note: in baseline test- both stay & switch pools have the same current
        %         %this is not true for the normal switching simulation..
        %
        
        
    case 'JK'
        
        %---baseline code---------
        
        
        if current_info.timeidx <= current_info.initpulse_time %initial pulse window, force first switch
            %stay cells get normal base current
            Iapp(Estay) = Iapp(Estay) + current_info.basecurr;
            %switch cells get half current, force first switch
            Iapp(Eswitch) = Iapp(Eswitch) + (.5 * current_info.basecurr);
        else %past initial pulse window
            if mod(numel(durations{state.stay}),2) == 0 %stim A (count is one behind, based off already recorded duration..)
                %all excitables get the same current
                Iapp(celltype.excit) = Iapp(celltype.excit) + (current_info.basecurr * current_info.stimA);
            else %stim B
                %all excitables get the same current
                Iapp(celltype.excit) = Iapp(celltype.excit) + (current_info.basecurr * current_info.stimB);
                %both excitatory & inhibitory cells get same current in JK model
            end
        end
        
        
        %note: in baseline test- both stay & switch pools have the same current
        %this is not true for the normal switching simulation..
        
        
end



