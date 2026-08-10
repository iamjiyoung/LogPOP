% Screen the candidate datasets used for Example 12(iv).
%
% J. Choi, July 31, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

candidate_seeds = (20261001:20261030)';
num_candidates = numel(candidate_seeds);
sample_size = 500;
num_starts = 20;
em_options.max_iterations = 1000;
em_options.relative_tolerance = 1e-7;

seed = string(candidate_seeds);
best_fraction = NaN(num_candidates, 1);
best_loglik = NaN(num_candidates, 1);
second_mode_gap = NaN(num_candidates, 1);
counts = strings(num_candidates, 1);

for i = 1:num_candidates
    fprintf('\n[%d/%d] Screening dataset seed %d\n', ...
        i, num_candidates, candidate_seeds(i));

    data = Ex_12_iv_generate_dataset(sample_size, candidate_seeds(i));
    em = Ex_12_iv_fit_lcm_em(data.patterns, data.counts, ...
        num_starts, candidate_seeds(i)+100000, em_options);

    rounded_values = unique(round(em.final_loglik, 5), 'sorted');
    num_best = sum(abs(em.final_loglik-em.loglik) <= 1e-5);
    best_fraction(i) = num_best/num_starts;
    best_loglik(i) = em.loglik;
    if numel(rounded_values) >= 2
        second_mode_gap(i) = rounded_values(end)-rounded_values(end-1);
    end
    counts(i) = string(mat2str(data.counts'));
end

results = table(seed, best_fraction, best_loglik, ...
    second_mode_gap, counts);
results = sortrows(results, {'best_fraction','second_mode_gap'}, ...
    {'ascend','descend'});

fprintf('\n=== Example 12(iv): dataset-screening results ===\n');
disp(results)
fprintf('Selected stress-instance seed : %s\n', results.seed(1));
