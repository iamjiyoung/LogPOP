function  comb = binom(N,m)

%  retun the number of combinations of choosing m out of N things 

if m == 0
   comb = 1;
   return;
end

comb = N/m;
for k=m-1:-1:1;
 comb = comb * (N+k-m)/k;
end
comb = round(comb);

return;