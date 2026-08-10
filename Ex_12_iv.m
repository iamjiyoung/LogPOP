% non-saturated latent-class-model results in Example 12(iv).
%
% J. Choi, August 10, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

data = Ex_12_iv_generate_dataset(500, 20261004);
em_options.max_iterations = 5000;
em_options.relative_tolerance = 1e-10;
fprintf('\n============================================================\n');
fprintf('Example 12(iv): 100-start EM calculation\n');
fprintf('============================================================\n');
em = Ex_12_iv_fit_lcm_em(data.patterns, data.counts, 100, 20361004, em_options);
fprintf('\n============================================================\n');
fprintf('Example 12(iv): order-3 moment relaxation\n');
fprintf('============================================================\n');
moment = Ex_12_iv_solve_moment(data.counts);

positive = data.counts > 0;
saturated_loglik = sum(data.counts(positive).* ...
    log(data.counts(positive)/sum(data.counts)));
moment_minus_saturated = moment.bound-saturated_loglik;
model_gap = saturated_loglik-em.loglik;
moment_gap = moment.bound-em.loglik;
deviance = 2*model_gap;
em_gaps = em.loglik-em.final_loglik;
num_max_iterations = sum(em.iterations == em_options.max_iterations);
largest_shortfall = max(em_gaps);

fprintf('\n=== Example 12(iv): reported results ===\n');
fprintf('Observed-variable count d       : %d\n', 4);
fprintf('Number of latent classes T      : %d\n', 2);
fprintf('Sample size                      : %d\n', sum(data.counts));
fprintf('Outcome counts                   : %s\n', mat2str(data.counts'));
fprintf('Saturated log-likelihood         : %.10f\n', saturated_loglik);
fprintf('Moment upper bound               : %.10f\n', moment.bound);
fprintf('Moment bound minus saturated     : %.3e\n', moment_minus_saturated);
fprintf('Best 100-start EM log-likelihood : %.10f\n', em.loglik);
fprintf('Moment bound minus EM            : %.10f\n', moment_gap);
fprintf('Deviance of best EM fit           : %.10f\n', deviance);
fprintf('Moment-relaxation runtime (seconds): %.4f\n', moment.time);
fprintf('100-start EM runtime (seconds)    : %.4f\n', em.time);
fprintf('Distinct EM terminal values      : %d\n', em.num_distinct_solutions);
fprintf('Starts within 1e-6 of best EM    : %d\n', em.num_best_solutions);
fprintf('Starts reaching 5000 iterations  : %d\n', num_max_iterations);
fprintf('Largest shortfall from best EM    : %.10f\n', largest_shortfall);
