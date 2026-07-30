function out = Ex_11_solve_case(report_weights, P, ord)
%EX_11_SOLVE_CASE Solve one paternity-analysis instance from Example 11.
% This helper is called repeatedly by Ex_11.m. Do not run it directly;
% Ex_11.m defines the data and reports the manuscript quantities. The
% output contains the LME bound, runtime, flat-truncation data, and atoms.
mset clear
mset('verbose', true)
nred = size(P,2)-1;
mpol('x', nred)

for i = 1:size(P,1)
    p{i} = P(i,end); %#ok<AGROW>
    for j = 1:nred
        p{i} = p{i} + (P(i,j)-P(i,end))*x(j);
    end
end

solve_weights = report_weights/max(report_weights);
for i = 1:numel(p)
    product = 1;
    for j = 1:numel(p)
        if j ~= i
            product = product*p{j};
        end
    end
    others{i} = product; %#ok<AGROW>
end

for i = 1:nred
    g{i} = 0; %#ok<AGROW>
    for j = 1:numel(p)
        g{i} = g{i} + solve_weights(j)*others{j}*diff(p{j},x(i));
    end
end
lme0 = 0;
for i = 1:nred
    lme0 = lme0 + x(i)*g{i};
end
for i = 1:nred
    lme(i) = lme0-g{i};
end

Ktheta = [sum(x)-1 <= 0; lme0 >= 0; lme0*(1-sum(x)) == 0];
for i = 1:nred
    Ktheta = [Ktheta; x(i) >= 0; lme(i) >= 0; lme(i)*x(i) == 0]; %#ok<AGROW>
end

flat_degree = ceil((numel(p)+1)/2);
out = raw_log_moment_relaxation(x, p, solve_weights, report_weights, ...
    Ktheta, ord, ord, flat_degree, flat_degree);
end

function out = raw_log_moment_relaxation( ...
        x, p, solve_weights, report_weights, Ktheta, ord, ...
        flat_max, flat_min, flat_degree)
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
    sdpsettings('solver', 'mosek', 'verbose', 1));
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
        1e-4, flat_min, false, flat_degree);
end
end
