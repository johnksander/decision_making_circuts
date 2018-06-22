function [CI,H] = boot_mean_diffs(x1,x2,Nstraps)
%boostrap test of mean differences 
alpha = .05;

x1 = bootstrp(Nstraps,@mean,x1);
x2 = bootstrp(Nstraps,@mean,x2);
strapped_diffs = x1 - x2;

CI = prctile(strapped_diffs,[100*alpha/2,100*(1-alpha/2)]);
if CI(1)>0 | CI(2)<0 %hypothesis test 
    H = 'reject H0';
else 
    H = 'fail to reject H0';
end


