function out = Ex_6_ii_solve_case(weights, use_lme, ord)
%EX_6_II_SOLVE_CASE Solve one coefficient instance from Example 6(ii).
% This helper is called by Ex_6_ii.m for either the standard or the LME
% moment relaxation. Do not run it directly; open and run Ex_6_ii.m.
% The output contains the relaxation value, runtime, flat-truncation data,
% extracted atoms, and the observed atom-to-bound gap.
%
% J. Choi, July 29, 2026
mset clear
mset('verbose', true)
mpol('x', 3)

p{1} = x(1)^3 + 3*x(1)^2*x(2) + 3*x(1)^2*x(3);
p{2} = 3*x(1)*x(2)^2 + 6*x(1)*x(2)*x(3);
p{3} = 3*x(1)*x(3)^2;
p{4} = x(2)^3 + 3*x(2)^2*x(3);
p{5} = 3*x(2)*x(3)^2;
p{6} = x(3)^3;

log_floor = 1e-5;
Ktheta = [x'*x <= 1; 1-sum(x) >= 0];
for i = 1:length(x)
    Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
end
for i = 1:length(p)
    Ktheta = [Ktheta; p{i} >= log_floor]; %#ok<AGROW>
end

flat_degree = 2;
if use_lme
    flat_degree = 3;
    g{1} = weights(1)*3*(x(1)+2*x(2)+2*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) ...
        +(weights(2)+weights(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));
    g{2} = weights(1)*3*x(2)*(x(2)+2*x(3))*(x(2)+3*x(3)) ...
        +weights(2)*(2*x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3)) ...
        +weights(4)*3*(x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3)) ...
        +weights(5)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));
    g{3} = weights(1)*3*x(3)*(x(2)+2*x(3))*(x(2)+3*x(3)) ...
        +weights(2)*2*x(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3)) ...
        +2*weights(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) ...
        +weights(4)*3*x(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3)) ...
        +(2*weights(5)+3*weights(6))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));
    lme0 = g{1}+g{2}+g{3};
    lme(1) = lme0*x(1)-g{1};
    lme(2) = lme0*x(2)-g{2};
    lme(3) = lme0*x(3)-g{3};
    Ktheta = [Ktheta; lme0 >= 0; lme0*(1-sum(x)) == 0];
    for i = 1:length(x)
        Ktheta = [Ktheta; lme(i) >= 0; lme(i)*x(i) == 0]; %#ok<AGROW>
    end
end

out = raw_log_moment_relaxation(x, p, weights, weights, Ktheta, ...
    ord, ord, flat_degree, flat_degree);
out.gap = NaN;
out.atom_values = NaN(1, numel(out.atoms));
for i = 1:numel(out.atoms)
    value_i = 0;
    for j = 1:numel(p)
        value_i = value_i + weights(j)*log(double(subs(p{j}, x, out.atoms{i})));
    end
    out.atom_values(i) = value_i;
end
if ~isempty(out.atom_values)
    out.gap = min(out.bound-out.atom_values);
end
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
