function  m = mon_order(power);

% input: 
%   pow: a interger vector of momomial power
% output: 
%      m: the graded lexicorgraphical ordering of 
%      monomial with power 
%      in all the monomials of degree <= d.

n = length(power);
d = sum(power);
if d > 0
m = binom(n+d-1,d-1) + hom_mon_order(power);
elseif d == 0
m = hom_mon_order(power);
end

