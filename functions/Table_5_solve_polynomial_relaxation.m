function out = Table_5_solve_polynomial_relaxation(x, objective_poly, Ktheta, ord, flat_max)
%TABLE_5_SOLVE_POLYNOMIAL_RELAXATION Solve one polynomial Moment-SOS relaxation.
% This lower-level helper is called by Table_5_solve_direct_product.m,
% not by a user-facing script. It returns the polynomial bound, runtime,
% moment-matrix ranks, flat-truncation order, and extracted atoms.
% Solve a polynomial moment relaxation with a linear moment objective.
%
% J. Choi, July 29, 2026

if nargin < 5
    flat_max = 0;
end

problem = msdp(max(objective_poly), Ktheta, ord);
[A, bp, c, cone, b0] = msedumi(problem);
mom_y = sdpvar(length(bp), 1);
c_Aty = c-A'*mom_y;
constraints = [];

for i = 1:cone.f
    constraints = [constraints, c_Aty(i) == 0]; %#ok<AGROW>
end
offset = cone.f;
for i = 1:cone.l
    constraints = [constraints, c_Aty(offset+i) >= 0]; %#ok<AGROW>
end
offset = offset+cone.l;

moment_sdp = [];
for i = 1:numel(cone.s)
    block_size = cone.s(i)^2;
    locM = mat(c_Aty(offset+1:offset+block_size));
    if i == 1
        moment_sdp = locM;
    end
    constraints = [constraints, locM >= 0]; %#ok<AGROW>
    offset = offset+block_size;
end

settings = sdpsettings('solver', 'mosek', 'verbose', 0);
objective = -(b0+bp'*mom_y);
warmup = optimize(constraints, objective, settings);
if warmup.problem ~= 0
    error('MOSEK warm-up failed: %s', warmup.info);
end
started = tic;
sol = optimize(constraints, objective, settings);
out.time = toc(started);
if sol.problem ~= 0
    error('MOSEK failed: %s', sol.info);
end

out.order = ord;
out.bound = value(b0+bp'*mom_y);
out.moment = value(moment_sdp);
out.full_rank = rank(out.moment, 1e-4);
out.flat_order = NaN;
out.flat_rank = NaN;
out.atoms = {};
if flat_max > 0
    previous_rng = rng;
    restore_rng = onCleanup(@() rng(previous_rng));
    rng(0, 'twister');
    [out.atoms, out.flat_order, out.flat_rank] = ...
        extractmin_flat(out.moment, length(x), flat_max, 1e-4, 1, false);
end
end
