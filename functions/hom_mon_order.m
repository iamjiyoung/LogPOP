function  m = hom_mon_order(power);

% input: 
%   power: a interger vector of momomial power
% output: 
%      m: the ordering of the poweromial  power 
%         in all the poweromials of degree d.

if length(find(power < 0)) > 0
    error(' the input vector of poweromial power must be nonnegative ');
end

n = length(power);
d = sum(power);

if n == 1 | d==0
   m = 1; 
   return;
end

m = 0;

nzidx = min(find(power));

if  nzidx >= 2
    zn = nzidx - 1;
    for k = d:-1:1
      m = m + binom(zn+k-1,k)*binom(n-zn+d-k-1,d-k);
    end  
      m = m + hom_mon_order(power(nzidx:n));
elseif nzidx == 1
    for k = d:-1:power(1)+1
      m = m + binom(n-1+d-k-1,d-k);
    end  
      m = m + hom_mon_order(power(2:n));
end

%
%  m = 1;
%
% if  power(1) == 0 
%     for k = d-1:-1:1
%     m = m + binom(n-1+d-k-1,d-k);
%     %m = m + nchoosek(n-1+d-k-1,d-k);
%     end   
%     if sum(power(2:n)) > 0
%      m = m + hom_mon_order(power(2:n));
%     % m = m + morder(power(2:n));
%     end
% elseif power(1) > 0
%     for k = d-1:-1:power(1)+1
%     m = m + binom(n-1+d-k-1,d-k);
%     %m = m + nchoosek(n-1+d-k-1,d-k);
%     end  
%     if sum(power(2:n)) > 0
%      % m = m + morder(power(2:n));
%      m = m + hom_mon_order(power(2:n));
%     end
% end

 