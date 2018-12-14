function lagged_var = next_timepoint(input_var)
%takes two column equation variable (T-1,T+0)
%sets T-1 = T+0, blanks T+0 for next timepoint. 
%"time-lagging" the variable. Fixed for 3d array input

sz = size(input_var);
lagged_var = NaN(sz);

if numel(sz) == 2
    lagged_var(:,1) = input_var(:,2);
elseif numel(sz) == 3
    lagged_var(:,1,:) = input_var(:,2,:);
end
