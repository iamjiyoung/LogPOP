% latent-class-model results in Example 12(i)-(iii).
%
% J. Choi, July 31, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

true_probabilities = {
    [0.60 0.40], ...
    [0.46 0.54], ...
    [0.26 0.34 0.09 0.31]};

% Cases (i) and (ii) retain the fixed datasets used in the paper.
% Case (iii) is generated from its stated data-generating parameters.
previous_rng = rng;
rng(20261003, 'twister');
draws = rand(500, 1);
counts3 = histcounts(draws, [0 cumsum(true_probabilities{3})]);
rng(previous_rng);
counts = {[303 197], [227 273], counts3};
d_values = [1; 1; 2];
class_values = [2; 3; 2];
instance = (1:3)';
k = 2*ones(3,1);
bound = NaN(3,1);
time = NaN(3,1);
true_loglik = NaN(3,1);
N = strings(3,1);
first_moments = strings(3,1);
case_names = ["i","ii","iii"];

fprintf('\nUnreported warm-up solve: Example 12(i)\n');
Ex_12_i_to_iii_solve_case(1, counts{1});

for i = 1:3
    fprintf('\n============================================================\n');
    fprintf('Example 12(%s), LCM moment relaxation, order k=2\n', ...
        case_names(i));
    fprintf('============================================================\n');
    out = Ex_12_i_to_iii_solve_case(i, counts{i});
    bound(i) = out.bound;
    time(i) = out.time;
    true_loglik(i) = sum(counts{i} .* log(true_probabilities{i}));
    N(i) = string(mat2str(counts{i}));
    first_moments(i) = string(mat2str(out.first_moments, 10));
end

fprintf('\n=== Example 12(i)-(iii): reported results ===\n');
for i = 1:3
    fprintf('\nCase (%s)\n', case_names(i));
    fprintf('Number of observed variables d : %d\n', d_values(i));
    fprintf('Number of latent classes T     : %d\n', class_values(i));
    fprintf('Outcome counts                 : %s\n', N(i));
    fprintf('Relaxation order k             : %d\n', k(i));
    fprintf('Moment upper bound             : %.10f\n', bound(i));
    fprintf('Log-likelihood at true params  : %.10f\n', true_loglik(i));
    fprintf('First moments                  : %s\n', first_moments(i));
    fprintf('Runtime (seconds)              : %.4f\n', time(i));
end
