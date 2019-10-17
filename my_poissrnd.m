function out = my_poissrnd(lambda,varargin)
%POISSRND Random arrays from the Poisson distribution.

% nope = lambda >= 15 | lambda < 0 | isinf(lambda);
% if any(nope),error('cannot handle lambda');end

out = arrayfun(@PFUNC, lambda);

function t = PFUNC(lambda)
t = 0;
p = 0 - log(rand());
while p < lambda
    t = t + 1; 
    p = p - log(rand());
end
