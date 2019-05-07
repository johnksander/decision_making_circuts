function val = get_percentile(X,perc)
%this gives better values than the matlab one, if I remember correctly
%percentile input should be in decimal format, ".1" gives top 10th percentile 

N = numel(X);
tind = [1:N] ./ N;
[~,tind] = min(abs(tind - (1 - perc)));
tind = min(tind); %in case there's a tie or something
X = sort(X); %watch out, this doubles X in memory 
val = X(tind);
