% paternity-analysis results in Example 11.
%
% J. Choi, August 10, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

Ndata = {
    [77 23], [63 37], [49 40 11], [83 2 15], [63 17 20], ...
    [59 8 16 17], [7 9 4 33 47], [39 38 23], [29 21 88 62]};
Pdata = {
    [0.5 1; 0.5 0], ...
    [0.5 1; 0.5 0], ...
    [0.5 0.875; 0.25 0.125; 0.25 0], ...
    [0.5 0.25; 0.5 0.5; 0 0.25], ...
    [0.5 0.25; 0.5 0.5; 0 0.25], ...
    [0.25 0.5; 0.25 0.5; 0.25 0; 0.25 0], ...
    [0.25 0; 0.25 0; 0.25 0; 0.25 0.5; 0 0.5], ...
    [0.5 0 0.875; 0.25 0.75 0.125; 0.25 0.25 0], ...
    [0.5 0 0; 0.25 0.75 0.25; 0.25 0.25 0.5; 0 0 0.25]};
orders = [2; 2; 3; 2; 2; 3; 4; 3; 4];

fprintf('\n============================================================\n');
fprintf('Example 11: unreported standard and LME warm-up solves\n');
fprintf('============================================================\n');
Ex_11_solve_case(Ndata{1}, Pdata{1}, 1, false);
Ex_11_solve_case(Ndata{1}, Pdata{1}, orders(1), true);

instance = (1:9)';
n = zeros(9,1);
g = zeros(9,1);
k = orders;
r = NaN(9,1);
t = NaN(9,1);
bound = NaN(9,1);
time = NaN(9,1);
standard_bound = NaN(9,1);
standard_time = NaN(9,1);
standard_optimizer = strings(9,1);
N = strings(9,1);
P = strings(9,1);
optimizer = strings(9,1);

for i = 1:9
    fprintf('\n============================================================\n');
    fprintf('Example 11, paternity instance %d of 9\n', i);
    fprintf('Standard moment relaxation, order k=1\n');
    fprintf('============================================================\n');
    standard = Ex_11_solve_case(Ndata{i}, Pdata{i}, 1, false);
    standard_bound(i) = standard.bound;
    standard_time(i) = standard.time;
    reduced_standard = standard.first_moments;
    full_standard = [reduced_standard, 1-sum(reduced_standard)];
    standard_optimizer(i) = string(mat2str(full_standard, 10));

    fprintf('\n');
    fprintf('LME moment relaxation, order k=%d\n', orders(i));
    fprintf('============================================================\n');
    out = Ex_11_solve_case(Ndata{i}, Pdata{i}, orders(i), true);
    n(i) = size(Pdata{i},2);
    g(i) = size(Pdata{i},1);
    r(i) = out.flat_rank;
    t(i) = out.flat_order;
    bound(i) = out.bound;
    time(i) = out.time;
    N(i) = string(mat2str(Ndata{i}));
    P(i) = string(mat2str(Pdata{i}));
    if ~isempty(out.atoms)
        reduced = out.atoms{1}(:)';
        full_x = [reduced, 1-sum(reduced)];
        optimizer(i) = string(mat2str(full_x, 10));
    end
end

fprintf('\n=== Example 11: reported results ===\n');
for i = 1:9
    fprintf('\nInstance %d\n', instance(i));
    fprintf('Number of fathers n       : %d\n', n(i));
    fprintf('Number of genotypes g     : %d\n', g(i));
    fprintf('Counts N                  : %s\n', N(i));
    fprintf('Standard order-1 value    : %.10f\n', standard_bound(i));
    fprintf('Standard first moments    : %s\n', standard_optimizer(i));
    fprintf('Standard runtime (seconds): %.4f\n', standard_time(i));
    fprintf('Relaxation order k        : %d\n', k(i));
    fprintf('LME relaxation value      : %.10f\n', bound(i));
    fprintf('Flat-truncation order t   : %d\n', t(i));
    fprintf('Rank                       : %d\n', r(i));
    fprintf('Extracted optimizer        : %s\n', optimizer(i));
    fprintf('Runtime (seconds)          : %.4f\n', time(i));
end

fprintf('\nMaximum absolute difference between standard and LME values: %.3e\n', ...
    max(abs(standard_bound-bound)));
