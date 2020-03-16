function [state,avail_noise] = check_noise_avail(stim_info,state)
%check and see what noise input the animal has available

avail_noise.Estay = 1;
avail_noise.Eswitch = 1;

if state.last_leave_end == state.timeidx
    %we just exited a leave state
    
    switch stim_info.delivery
        case 'constant'
            
            state.sample_clock = 0;
            
        case 'pulse'
            
            
            
            
            switch stim_info.sample_schedule
                case 'flexible' %sample schedule based on behavior
                    
                    switch state.state_def
                        case 'include_undecided'
                            %this was intended for use with an undecided delay period!!
                            
                            %set the sample timing to -delay period
                            state.sample_clock = -stim_info.pulse(2); %this adds an extra delay here
                            
                        case 'active_states'
                            state.sample_clock = 0;
                            
                    end
                case 'fixed'
                    %do nothing, just stick to the schedule 
            end
    end
end



%advance the sample clock
state.sample_clock = state.sample_clock + 1;



if all(state.switch == state.now)
    %we're in a leave state
    if state.count >= state.cut_leave_state
        %it's been X ms, we've left and now we're kicking off the leave state
        
        switch state.state_def  %whether simulation aknowledges "undecided states"
            case 'active_states'
                %force a transition to stay
                avail_noise.Eswitch = .5;
            case 'include_undecided'
                %force an undecided, half noise to all E cells
                avail_noise.Estay = .5;
                avail_noise.Eswitch = .5;
        end
    end
    
    
else %we're either in a stay, or undecided state.
    
    
    switch stim_info.delivery
        case 'constant'
            %nothing special to do here
            
        case 'pulse'
            
            switch state.state_def
                case 'include_undecided' %below code shouldn't really be used if you're not including undecided
                    
                    %Determine what part of the sampling pulse we're in:
                    %----If we're in the off-part of a pulse, half to all E-cell (undecided)
                    %----If the sample clock is negative, it's the delay after leave (also force undecided)
                    %Sample_clock was set to -delay after the last "transition out of leave" (to undecided, presumably).
                    %----If the sample clock is positive, and < E-stay_offset, half to E-switch
                    %kick into the stay state. You only do this after the bistability check passes
                    
                    seq_length = sum(stim_info.pulse); %how long for a single on, off sequence
                    t = mod(state.sample_clock,seq_length); %find if we're in the on, or off part of that sequence
                    if t > stim_info.pulse(1) || state.sample_clock < 1
                        %we're past the "stimulus on" part, the stimulus is off
                        %or, the sample_clock is negative: we're in the delay after leave
                        avail_noise.Estay = .5;
                        avail_noise.Eswitch = .5;
                    elseif t <= state.sample_Estay_offset && state.sample_clock > 0
                        %beginning part of the sample, kick on stay-state
                        avail_noise.Eswitch = .5;
                    end
            end
    end
end