function instance = Table_5_6_generate_instance(n, m, seed, epsilon)
%TABLE_5_6_GENERATE_INSTANCE Generate a fixed quadratic test instance.
% This helper is called by Table_5.m and Table_6.m. Do not run it
% directly; the table scripts set the dimensions, seed, and epsilon.
% The returned structure stores all coefficients and logarithmic weights.
% Generate a reproducible positive, nonconcave quadratic log-POP instance.
%
% J. Choi, July 29, 2026

if nargin < 4
    epsilon = 0.05;
end
if n < 1 || m < 1 || fix(n) ~= n || fix(m) ~= m
    error('n and m must be positive integers.');
end
if epsilon <= 0
    error('epsilon must be positive.');
end

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(seed, 'twister');

U = zeros(m, n);
beta = zeros(m, 1);
for j = 1:m
    U(j, :) = randn(1, n);
    U(j, :) = U(j, :)/norm(U(j, :));
    beta(j) = 0.15*(2*rand-1);
end

instance.n = n;
instance.m = m;
instance.seed = seed;
instance.epsilon = epsilon;
instance.U = U;
instance.beta = beta;
instance.weights = ones(m, 1);
end
