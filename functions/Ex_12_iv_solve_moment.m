function out = Ex_12_iv_solve_moment(report_weights)
%EXAMPLE_12_IV_SOLVE_MOMENT Solve the moment relaxation for Example 12(iv).
% This helper is called by Ex_12_iv.m and should not be run directly.
% The output contains the moment upper bound, runtime, moment matrix,
% first moments, and flat-truncation information when available.
%
% J. Choi, August 10, 2026
mset clear
mset('verbose', true)
mpol('x', 9)

patterns = dec2bin(0:15)-'0';
p = cell(16, 1);
for r = 1:16
    component1 = 1;
    component2 = 1;
    for j = 1:4
        component1 = component1*x(1+j)^patterns(r,j)* ...
            (1-x(1+j))^(1-patterns(r,j));
        component2 = component2*x(5+j)^patterns(r,j)* ...
            (1-x(5+j))^(1-patterns(r,j));
    end
    p{r} = x(1)*component1+(1-x(1))*component2;
end

Ktheta = [9-x'*x >= 0];
for i = 1:9
    Ktheta = [Ktheta; x(i) >= 0; 1-x(i) >= 0]; %#ok<AGROW>
end
solve_weights = report_weights/max(report_weights);
out = raw_log_moment_relaxation(x, p, solve_weights, report_weights, ...
    Ktheta, 3, 0, 1);
end

function out = raw_log_moment_relaxation( ...
        x, p, solve_weights, report_weights, Ktheta, ord, flat_max, flat_min)
% Build and solve the moment relaxation explicitly.

% Step 1: Build a GloptiPoly moment relaxation for each p_i.
num_terms = numel(p);
PO = cell(num_terms, 1);
A = cell(num_terms, 1);
bp = cell(num_terms, 1);
c = cell(num_terms, 1);
K = cell(num_terms, 1);
b0 = cell(num_terms, 1);
for i = 1:num_terms
    PO{i} = msdp(max(p{i}), Ktheta, ord);
    [A{i}, bp{i}, c{i}, K{i}, b0{i}] = msedumi(PO{i});
end

% Step 2: Introduce the truncated moment variables.
moment_cone = K{1};
y = sdpvar(length(bp{1}), 1);
cone_coordinates = c{1}-A{1}'*y;
constraints = [];

% Step 3: Recover the equality, inequality, and PSD constraints.
for i = 1:moment_cone.f
    constraints = [constraints, cone_coordinates(i) == 0]; %#ok<AGROW>
end

offset = moment_cone.f;
for i = 1:moment_cone.l
    constraints = [constraints, cone_coordinates(offset+i) >= 0]; %#ok<AGROW>
end
offset = offset+moment_cone.l;

moment_matrix = [];
for i = 1:numel(moment_cone.s)
    block_size = moment_cone.s(i)^2;
    psd_block = mat(cone_coordinates(offset+1:offset+block_size));
    constraints = [constraints, psd_block >= 0]; %#ok<AGROW>
    if i == 1
        moment_matrix = psd_block;
    end
    offset = offset+block_size;
end

% Step 4: Form the logarithmic objective.
log_arguments = cell(num_terms, 1);
objective = 0;
for i = 1:num_terms
    log_arguments{i} = b0{i}+bp{i}'*y;
    objective = objective-solve_weights(i)*log(log_arguments{i});
end

% Step 5: Solve the SDP with MOSEK.
started = tic;
solution = optimize(constraints, objective, ...
    sdpsettings('solver', 'mosek', 'verbose', 0));
elapsed = toc(started);
if solution.problem ~= 0
    error('MOSEK failed: %s', solution.info);
end

% Step 6: Read the bound and the optimal moment matrix.
arguments = zeros(1, num_terms);
for i = 1:num_terms
    arguments(i) = value(log_arguments{i});
end
out.order = ord;
out.bound = sum(report_weights(:)'.*log(arguments));
out.time = elapsed;
out.arguments = arguments;
out.first_moments = value(y(1:length(x)))';
out.moment = value(moment_matrix);
out.full_rank = rank(out.moment, 1e-4);
out.flat_order = NaN;
out.flat_rank = NaN;
out.atoms = {};

% Step 7: Check flat truncation and extract atoms when requested.
if flat_max > 0
    previous_rng = rng;
    restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
    rng(0, 'twister');
    [out.atoms, out.flat_order, out.flat_rank] = ...
        extractmin_flat(out.moment, length(x), flat_max, ...
        1e-4, flat_min, false);
end
end
