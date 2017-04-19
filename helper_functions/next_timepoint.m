function lagged_var = next_timepoint(input_var)
%takes two column equation variable (T-1,T+0)
%sets T-1 = T+0, blanks T+0 for next timepoint. 
%"time-lagging" the variable. 

lagged_var = NaN(size(input_var));
lagged_var(:,1) = input_var(:,2);
