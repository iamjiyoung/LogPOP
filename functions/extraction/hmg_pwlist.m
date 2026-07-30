function  pwlst = hmg_pwlist(nvar,deg);

% input: 
%   nvar: numbers of variables
%    deg: degree of the monomials
% output: 
%  pwlst: an integer matrix whose rows are monomials powers 
%         of the same degree deg. 
%        The rows are ordered from top to bottom lexicographically
%
 

if deg == 0 
    pwlst = zeros(1,nvar);
    return;
elseif deg == 1
    pwlst = eye(nvar);
    return;
elseif nvar == 1
    pwlst = [deg];  
    return;
end

pwlst = [deg  zeros(1,nvar-1)];

for k = deg-1: -1 : 0
   A = hmg_pwlist(nvar-1,deg-k);
   pwlst = [ pwlst; k*ones(size(A,1),1) A];
end

return;