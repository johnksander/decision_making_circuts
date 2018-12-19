function segIDs = get_overlapIDs(segIDs,dup_inds)
%Set all overlapping segments to the same seqID. Diff will show the continously overlapping seq IDs.
%There may be overlap > 2 IDs, so you'll want to do the full overlapping sequence!


overlap_seq = segIDs(dup_inds);

%this distinguishes different sets of overlapping segments 
seg_breaks = [0;diff(overlap_seq)>1];

curr_seg = overlap_seq(1); %start with the first one
for idx = 1:numel(overlap_seq)
    %check and see if we're in a new segment or not
    if seg_breaks(idx) == 1
        curr_seg = overlap_seq(idx); %pick up the new starting ID
    end
    fixIDs = overlap_seq(idx);
    segIDs(segIDs == fixIDs) = curr_seg; %assign overlapped segIDs to the first ID in that segment 
end

