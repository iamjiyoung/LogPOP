function M=mon_pow_vec(n,d)

% M is the matrix whose rows are powers of 
% ALL monomials of n indeterminates and degree <= d
% ordered graded lexicographically

M = [];

for deg = 0: 1 : d
  M = [M; hmg_pwlist(n,deg)];
end
