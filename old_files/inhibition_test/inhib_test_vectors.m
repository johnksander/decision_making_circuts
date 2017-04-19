function [Eon,Ioff] = inhib_test_vectors(results)

muEstay_Sg = results.muEstay_Sg;
muEswitch_Sg = results.muEswitch_Sg;
Erate_stay = results.Erate_stay;
Erate_switch = results.Erate_switch;
Irate_stay = results.Irate_stay;
Irate_switch = results.Irate_switch;

sz = size(Irate_switch);
num_timepoints = numel(Irate_switch);
Eon = NaN(sz);
Ioff = NaN(sz);


Irates = [Irate_stay;Irate_switch];
Erates = [Erate_stay;Erate_switch];
Sg = [muEstay_Sg;muEswitch_Sg];


start_idx = muEstay_Sg *4 < muEswitch_Sg | muEswitch_Sg*4 < muEstay_Sg;
start_idx = find(start_idx,1,'first');
curr_state = Sg(:,start_idx) == max(Sg(:,start_idx));

for idx = start_idx:num_timepoints
    if Sg(curr_state,idx) * 4 < Sg(~curr_state,idx) 
        %switch!
        curr_state = ~curr_state;
    end
    Eon(idx) = Erates(curr_state,idx);
    Ioff(idx) = Irates(~curr_state,idx);
end
