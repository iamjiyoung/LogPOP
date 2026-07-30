% latent-class-model results in Example 12(i)-(iii).
%
% J. Choi, July 29, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

counts = {[303 197], [227 273], [97 211 100 92]};
d_values = [1; 1; 2];
class_values = [2; 3; 2];
instance = (1:3)';
k = 2*ones(3,1);
bound = NaN(3,1);
time = NaN(3,1);
N = strings(3,1);
first_moments = strings(3,1);
case_names = ["i","ii","iii"];

for i = 1:3
    fprintf('\n============================================================\n');
    fprintf('Example 12(%s), LCM moment relaxation, order k=2\n', ...
        case_names(i));
    fprintf('============================================================\n');
    out = Ex_12_i_to_iii_solve_case(i, counts{i});
    bound(i) = out.bound;
    time(i) = out.time;
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
    fprintf('First moments                  : %s\n', first_moments(i));
    fprintf('Runtime (seconds)              : %.4f\n', time(i));
end
