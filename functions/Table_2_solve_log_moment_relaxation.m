function out = Table_2_solve_log_moment_relaxation( ...
        x, p, solve_weights, report_weights, Ktheta, ord, ...
        flat_max, flat_min, flat_degree)
% Solve one log-moment relaxation for Table 2.
%
% J. Choi, July 29, 2026

if nargin < 7
    flat_max = 0;
end
if nargin < 8
    flat_min = 1;
end
if nargin < 9
    flat_degree = 1;
end

solve_weights = solve_weights(:)';
report_weights = report_weights(:)';
if numel(p) ~= numel(solve_weights) || numel(p) ~= numel(report_weights)
    error('The number of weights must equal the number of log terms.');
end

for i = 1:numel(p)
    PO{i} = msdp(max(p{i}), Ktheta, ord); %#ok<AGROW>
    [A{i}, bp{i}, c{i}, K{i}, b0{i}] = msedumi(PO{i}); %#ok<AGROW>
end

MCone = K{1};
mom_y = sdpvar(length(bp{1}), 1);
c_Aty = c{1} - A{1}'*mom_y;
MomRelax = [];

for i = 1:MCone.f
    MomRelax = [MomRelax, c_Aty(i) == 0]; %#ok<AGROW>
end

offset = MCone.f;
for i = 1:MCone.l
    MomRelax = [MomRelax, c_Aty(offset+i) >= 0]; %#ok<AGROW>
end
offset = offset + MCone.l;

moment_sdp = [];
for i = 1:numel(MCone.s)
    block_size = MCone.s(i)^2;
    locM = mat(c_Aty(offset+1:offset+block_size));
    if i == 1
        moment_sdp = locM;
    end
    MomRelax = [MomRelax, locM >= 0]; %#ok<AGROW>
    offset = offset + block_size;
end

obj = 0;
for i = 1:numel(p)
    obj = obj - solve_weights(i)*log(b0{i} + bp{i}'*mom_y);
end

settings = sdpsettings('solver', 'mosek', 'verbose', 0);
warmup = optimize(MomRelax, obj, settings);
if warmup.problem ~= 0
    error('MOSEK warm-up failed: %s', warmup.info);
end
started = tic;
sol = optimize(MomRelax, obj, settings);
elapsed = toc(started);
if sol.problem ~= 0
    error('MOSEK failed: %s', sol.info);
end

arguments = zeros(1, numel(p));
for i = 1:numel(p)
    arguments(i) = value(b0{i} + bp{i}'*mom_y);
end

out.order = ord;
out.bound = sum(report_weights .* log(arguments));
out.time = elapsed;
out.arguments = arguments;
out.first_moments = value(mom_y(1:length(x)))';
out.moment = value(moment_sdp);
out.full_rank = rank(out.moment, 1e-4);
out.flat_order = NaN;
out.flat_rank = NaN;
out.atoms = {};

if flat_max > 0
    previous_rng = rng;
    restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
    rng(0, 'twister');
    [out.atoms, out.flat_order, out.flat_rank] = ...
        extractmin_flat(out.moment, length(x), flat_max, ...
        1e-4, flat_min, false, flat_degree);
end
end
