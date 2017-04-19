function Iapp = addpulse(pulse_params,Iapp,timevec)
%add current pulse to Iapp, specified in pulse params

num_pulses = numel(pulse_params(:,1));
for idx = 1:num_pulses
    
    pulse_targets = pulse_params{idx,1};
    pulse_current = pulse_params{idx,2};
    pulse_start = pulse_params{idx,3};
    pulse_end = pulse_params{idx,4};
    
    pulse_timecourse = timevec >= pulse_start & timevec <= pulse_end;
    Iapp(pulse_targets,pulse_timecourse) = Iapp(pulse_targets,pulse_timecourse) + pulse_current;
    
    
    
end