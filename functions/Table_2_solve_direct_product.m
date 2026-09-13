function out = Table_2_solve_direct_product(instance, min_ord, max_extra_orders, max_psd_dim)
% Solve the direct-product relaxation for Table 2.
%
% J. Choi, July 29, 2026

integer_weights = round(instance.weights);
if any(abs(instance.weights-integer_weights) > 1e-10) || any(integer_weights < 1)
    error('FullSOS:NonintegerWeights', ...
        'The direct polynomial formulation requires positive integer weights.');
end
required_order = sum(integer_weights);
if nargin < 2 || isempty(min_ord)
    min_ord = required_order;
end
if min_ord < required_order
    error('The full SOS order must be at least the sum of integer weights.');
end
if nargin < 3
    max_extra_orders = 2;
end
if nargin < 4
    max_psd_dim = Inf;
end

mset clear
mset('verbose', false)
mpol('x', instance.n)

objective_poly = 1;
for j = 1:instance.m
    affine = instance.beta(j);
    for i = 1:instance.n
        affine = affine+instance.U(j, i)*x(i);
    end
    objective_poly = objective_poly*(instance.epsilon+affine^2)^integer_weights(j);
end

Ktheta = [x'*x <= 1; 1-sum(x) >= 0];
for i = 1:instance.n
    Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
end

total_time = 0;
out = [];
for ord = min_ord:min_ord+max_extra_orders
    if nchoosek(instance.n+ord, ord) > max_psd_dim
        break
    end
    candidate = Table_2_solve_polynomial_relaxation(x, objective_poly, Ktheta, ord, ord);
    total_time = total_time+candidate.time;
    out = candidate;
    if ~isnan(candidate.flat_order)
        break
    end
end
if isempty(out)
    error('FullSOS:SizeLimit', 'No full SOS order satisfied the PSD-size limit.');
end
out.time = total_time;
end
