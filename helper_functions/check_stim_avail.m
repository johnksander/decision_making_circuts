function [state,outcome] = check_stim_avail(stim_info,state)
%check and see if the animal has a stimulus available
%----
%return false if the stimulus delivery method is "pulse"
%and we're currently in the "stimulus off" part of that delivery squence.
%(must be in the stay-state as well, for clarity).
%----
%return false if we've been in a leave state for Xms, (cut_leave_state)
%----
%return false if it's been less than ISI ("stimulus off" part) since the last leave
%state (don't check for last leave state, do this 'smartly' with a counter)
%----
%otherwise return true for any other scenario

%fix the description above dude 

outcome = true;

switch stim_info.delivery
    case 'pulse'
        
        if state.last_leave_end == state.timeidx
            %we just exited a leave state
            %set the sample timing to -delay period
            state.sample_clock = -stim_info.pulse(2); %this adds an extra delay here
        end
        
        %advance the sample clock
        state.sample_clock = state.sample_clock + 1;
        
        if all(state.switch == state.now)
            %we're in a leave state
            if state.count >= state.cut_leave_state
                %it's been X ms, we've left and the stimulus is no longer available
                outcome = false;
            end
        else
            %we're either in a stay, or undecided state. Sample_clock was
            %set to -delay after the last "transition out of leave" (to undecided,
            %presumably). Determine what part of the sampling pulse we're
            %in, if the sample clock is positive
            
            seq_length = sum(stim_info.pulse); %how long for a single on, off sequence
            t = mod(state.sample_clock,seq_length); %find if we're in the on, or off part of that sequence
            if t > stim_info.pulse(1) || state.sample_clock < 1
                %we're past the "stimulus on" part, the stimulus is off
                %Or, the sample_clock is negative: we're in the delay after leave 
                outcome = false;
            end
        end
end


%         if all(state.stay == state.now)
%             %we're in a stay state
%             seq_length = sum(stim_info.pulse); %how long for a single on, off sequence
%             t = mod(state.count,seq_length); %find if we're in the on, or off part of that sequence
%             if t > stim_info.pulse(1)
%                 %we're past the "stimulus on" part, the stimulus is off
%                 outcome = false;
%             end
%
%         elseif all(state.switch == state.now)
%             %we're in a leave state
%             if state.count >= state.cut_leave_state
%                 %it's been X ms, we've left and the stimulus is no longer available
%                 outcome = false;
%             end
%
%         elseif all(state.undecided == state.now)
%
%         elseif state.last_leave_end <= stim_info.pulse(2)
%             keyboard
%         end
