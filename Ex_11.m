% paternity-analysis results in Example 11.

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

instance = (1:9)';
n = zeros(9,1);
t = zeros(9,1);
k = orders;
r = NaN(9,1);
flat_t = NaN(9,1);
bound = NaN(9,1);
time = NaN(9,1);
N = strings(9,1);
P = strings(9,1);
optimizer = strings(9,1);

for i = 1:9
    fprintf('\n============================================================\n');
    fprintf('Example 11, paternity instance %d of 9\n', i);
    fprintf('LME moment relaxation, order k=%d\n', orders(i));
    fprintf('============================================================\n');
    out = Ex_11_solve_case(Ndata{i}, Pdata{i}, orders(i));
    n(i) = size(Pdata{i},2);
    t(i) = size(Pdata{i},1);
    r(i) = out.flat_rank;
    flat_t(i) = out.flat_order;
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
    fprintf('Number of genotypes t     : %d\n', t(i));
    fprintf('Counts N                  : %s\n', N(i));
    fprintf('Relaxation order k        : %d\n', k(i));
    fprintf('Optimal relaxation value  : %.10f\n', bound(i));
    fprintf('Flat-truncation order t   : %d\n', flat_t(i));
    fprintf('Rank                       : %d\n', r(i));
    fprintf('Extracted optimizer        : %s\n', optimizer(i));
    fprintf('Runtime (seconds)          : %.4f\n', time(i));
end
