% Reproduce the rank-two LME relaxation reported in Section 5.

clear
clc
mset clear
mpol x 1

p{1} = 2+x^2;
p{2} = 3+x^2;

% Box constraints and their polynomialized LME constraints.
Ktheta = [x+1 >= 0; 1-x >= 0];
g = 2*x*p{2}+2*x*p{1};
hat_tau_lower = g*(x-1)/2;
hat_tau_upper = g*(x+1)/2;
Ktheta = [Ktheta; ...
    hat_tau_lower >= 0; hat_tau_upper >= 0; ...
    hat_tau_lower*(x+1) == 0; ...
    hat_tau_upper*(1-x) == 0];

ord = 4;
for i = 1:length(p)
    problem{i} = msdp(max(p{i}), Ktheta, ord); %#ok<SAGROW>
    [A{i}, bp{i}, c{i}, K{i}, b0{i}] = ...
        msedumi(problem{i}); %#ok<SAGROW>
end

moment_cone = K{1};
moment_y = sdpvar(length(bp{1}), 1);
cone_coordinates = c{1}-A{1}'*moment_y;
constraints = [];

for i = 1:moment_cone.f
    constraints = [constraints, cone_coordinates(i) == 0]; %#ok<AGROW>
end

offset = moment_cone.f;
for i = 1:moment_cone.l
    constraints = [constraints, ...
        cone_coordinates(offset+i) >= 0]; %#ok<AGROW>
end
offset = offset+moment_cone.l;

moment_matrix = [];
for i = 1:length(moment_cone.s)
    block_size = moment_cone.s(i)^2;
    psd_block = mat(cone_coordinates(offset+1:offset+block_size));
    constraints = [constraints, psd_block >= 0]; %#ok<AGROW>
    if i == 1
        moment_matrix = psd_block;
    end
    offset = offset+block_size;
end

log_argument_1 = b0{1}+bp{1}'*moment_y;
log_argument_2 = b0{2}+bp{2}'*moment_y;
objective = -log(log_argument_1)-log(log_argument_2);

solution = optimize(constraints, objective, ...
    sdpsettings('solver', 'mosek', 'verbose', 1));
if solution.problem ~= 0
    error('MOSEK failed: %s', solution.info)
end

optimal_value = log(value(log_argument_1))+log(value(log_argument_2));
moment_matrix = value(moment_matrix);
rank_M4 = rank(moment_matrix, 1e-4);
rank_M1 = rank(moment_matrix(1:2, 1:2), 1e-4);
atoms = extractmin(moment_matrix, 1, ord);
extracted_points = sort(cellfun(@(point) point(1), atoms));

fprintf('\n=== Rank-two LME relaxation ===\n');
fprintf('Relaxation order k        : %d\n', ord);
fprintf('Optimal relaxation value  : %.10f\n', optimal_value);
fprintf('Difference from log(12)   : %.3e\n', optimal_value-log(12));
fprintf('Flat-truncation order t    : 4\n');
fprintf('Rank M_4                   : %d\n', rank_M4);
fprintf('Rank M_1                   : %d\n', rank_M1);
fprintf('Extracted maximizers       : %s\n', ...
    mat2str(extracted_points, 12));
