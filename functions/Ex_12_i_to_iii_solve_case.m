function out = Ex_12_i_to_iii_solve_case(instance, report_weights)
%EX_12_I_TO_III_SOLVE_CASE Solve one saturated LCM instance.
% This helper is called by Ex_12_i_to_iii.m for cases (i), (ii), and (iii).
% Do not run it directly. The output contains the moment upper bound,
% runtime, logarithm arguments, and first moments.
%
% J. Choi, July 29, 2026
mset clear
mset('verbose', true)

switch instance
    case 1
        mpol('x', 6)
        p{1} = x(1)*x(3) + x(2)*x(5);
        p{2} = x(1)*x(4) + x(2)*x(6);
        Ktheta = [x(1)+x(2) == 1; x(3)+x(4) == 1; x(5)+x(6) == 1];
    case 2
        mpol('x', 9)
        p{1} = x(1)*x(4) + x(2)*x(6) + x(3)*x(8);
        p{2} = x(1)*x(5) + x(2)*x(7) + x(3)*x(9);
        Ktheta = [x(1)+x(2)+x(3) == 1; x(4)+x(5) == 1; ...
            x(6)+x(7) == 1; x(8)+x(9) == 1];
    case 3
        mpol('x', 10)
        p{1} = x(1)*x(3)*x(5) + x(2)*x(7)*x(9);
        p{2} = x(1)*x(3)*x(6) + x(2)*x(7)*x(10);
        p{3} = x(1)*x(4)*x(5) + x(2)*x(8)*x(9);
        p{4} = x(1)*x(4)*x(6) + x(2)*x(8)*x(10);
        Ktheta = [x(1)+x(2) == 1; x(3)+x(4) == 1; ...
            x(5)+x(6) == 1; x(7)+x(8) == 1; x(9)+x(10) == 1];
    otherwise
        error('Unknown LCM instance.');
end

for i = 1:length(x)
    Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
end
solve_weights = report_weights/max(report_weights);
out = raw_log_moment_relaxation(x, p, solve_weights, report_weights, ...
    Ktheta, 2, 0, 1);
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
        1e-4, flat_min, false);
end
end
