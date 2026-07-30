% J. Choi, July 29, 2026
%
% non-saturated latent-class-model results in Example 12(iv).

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
sorted_em_gaps = sort(em_gaps);
em_p90_gap = sorted_em_gaps(ceil(0.9*numel(sorted_em_gaps)));

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
fprintf('Deviance from saturated model    : %.10f\n', deviance);
fprintf('Moment runtime (seconds)          : %.4f\n', moment.time);
fprintf('100-start EM runtime (seconds)    : %.4f\n', em.time);
fprintf('Distinct EM terminal values      : %d\n', em.num_distinct_solutions);
fprintf('Starts within 1e-6 of best EM    : %d\n', em.num_best_solutions);
