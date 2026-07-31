% the rank-two LME relaxation in Example 7.
%
% J. Choi, July 30, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

standard = [];
fprintf(['\nThe standard order-k=1 relaxation is unbounded for the ', ...
    'box-constraint representation used here.\n']);
for candidate_order = 2:3
    fprintf('\n--- Standard moment relaxation, order k=%d ---\n', ...
        candidate_order);
    try
        candidate = Ex_7_solve_case(false, candidate_order);
    catch solver_error
        fprintf('No optimal solution returned at k=%d: %s\n', ...
            candidate_order, solver_error.message);
        continue
    end
    if isnan(candidate.flat_order)
        fprintf('No flat truncation detected at k=%d (rank M_%d = %d).\n', ...
            candidate_order, candidate_order, candidate.full_rank);
    else
        standard = candidate;
        break
    end
end
if isempty(standard)
    error('No standard flat truncation found through order 3.');
end

lme = [];
for candidate_order = 3:4
    fprintf('\n--- LME moment relaxation, order k=%d ---\n', ...
        candidate_order);
    try
        candidate = Ex_7_solve_case(true, candidate_order);
    catch solver_error
        fprintf('No optimal solution returned at k=%d: %s\n', ...
            candidate_order, solver_error.message);
        continue
    end
    if isnan(candidate.flat_order)
        fprintf('No flat truncation detected at k=%d (rank M_%d = %d).\n', ...
            candidate_order, candidate_order, candidate.full_rank);
    else
        lme = candidate;
        break
    end
end
if isempty(lme)
    error('No LME flat truncation found through order 4.');
end
solutions = {standard; lme};

method = ["standard"; "LME"];
order = [standard.order; lme.order];
bound = [standard.bound; lme.bound];
true_value = log(12)*ones(2, 1);
bound_gap = bound-true_value;
solve_time = [standard.time; lme.time];
flat_order = [standard.flat_order; lme.flat_order];
flat_rank = [standard.flat_rank; lme.flat_rank];
atoms = strings(2, 1);
log_arguments = strings(2, 1);

rank_tol_1e_6 = NaN(2, 1);
rank_tol_1e_5 = NaN(2, 1);
rank_tol_1e_4 = NaN(2, 1);
rank_tol_1e_3 = NaN(2, 1);
tolerances = [1e-6, 1e-5, 1e-4, 1e-3];

for i = 1:2
    out = solutions{i};
    values = zeros(1, numel(out.atoms));
    for j = 1:numel(out.atoms)
        values(j) = out.atoms{j}(1);
    end
    atoms(i) = string(mat2str(sort(values), 12));
    log_arguments(i) = string(mat2str(out.arguments, 12));
    ranks = arrayfun(@(tol) rank(out.moment, tol), tolerances);
    rank_tol_1e_6(i) = ranks(1);
    rank_tol_1e_5(i) = ranks(2);
    rank_tol_1e_4(i) = ranks(3);
    rank_tol_1e_3(i) = ranks(4);
end

fprintf('\n=== Example 7: reported results ===\n');
for i = 1:2
    fprintf('\n%s moment relaxation\n', method(i));
    fprintf('Relaxation order k        : %d\n', order(i));
    fprintf('Optimal relaxation value  : %.10f\n', bound(i));
    fprintf('Difference from log(12)   : %.3e\n', bound_gap(i));
    fprintf('Flat-truncation order t    : %d\n', flat_order(i));
    fprintf('Rank                       : %d\n', flat_rank(i));
    fprintf('Extracted atoms            : %s\n', atoms(i));
    fprintf('Runtime (seconds)          : %.4f\n', solve_time(i));
end
fprintf('\nLME moment-matrix ranks by tolerance\n');
fprintf('Tolerance 1e-6            : %d\n', rank_tol_1e_6(2));
fprintf('Tolerance 1e-5            : %d\n', rank_tol_1e_5(2));
fprintf('Tolerance 1e-4            : %d\n', rank_tol_1e_4(2));
fprintf('Tolerance 1e-3            : %d\n', rank_tol_1e_3(2));
