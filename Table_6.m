% Reproduce the MATLAB results reported in Table 6.
% Run this script directly.
%
% J. Choi, July 29, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

config.n = 4;
config.m = 3;
config.seeds = 3;
config.epsilon = 0.05;
config.num_starts = 100;
config.profile_names = ["unit"; "small_integer"; "moderate_integer"; ...
    "large_integer"; "fractional"; "irrational"];
config.weight_profiles = [1 1 1; 1 2 3; 2 5 10; 5 10 20; ...
    0.5 1.5 2.5; sqrt(2) pi/2 exp(1)/2];
config.max_full_sos_psd_dim = 250;
config.full_sos_max_extra_orders = 2;

num_profiles = size(config.weight_profiles, 1);
num_cases = num_profiles*numel(config.seeds);
profile = strings(num_cases, 1);
seed = NaN(num_cases, 1);
n = config.n*ones(num_cases, 1);
m = config.m*ones(num_cases, 1);
weights_text = strings(num_cases, 1);
sum_weights = NaN(num_cases, 1);
integer_weights = false(num_cases, 1);
local_value = NaN(num_cases, 1);
local_time = NaN(num_cases, 1);
local_successes = NaN(num_cases, 1);
local_solver = strings(num_cases, 1);
mom_bound = NaN(num_cases, 1);
mom_time = NaN(num_cases, 1);
mom_psd_dim = (config.n+1)*ones(num_cases, 1);
mom_log_gap = NaN(num_cases, 1);
mom_gap_per_weight = NaN(num_cases, 1);
mom_multiplicative_gap = NaN(num_cases, 1);
lme_order = (config.m+1)*ones(num_cases, 1);
lme_psd_dim = nchoosek(config.n+config.m+1, config.m+1)*ones(num_cases, 1);
lme_bound = NaN(num_cases, 1);
lme_time = NaN(num_cases, 1);
lme_log_gap = NaN(num_cases, 1);
lme_gap_per_weight = NaN(num_cases, 1);
lme_flat_order = NaN(num_cases, 1);
lme_flat_rank = NaN(num_cases, 1);
full_sos_required_order = NaN(num_cases, 1);
full_sos_required_psd_dim = NaN(num_cases, 1);
full_sos_actual_order = NaN(num_cases, 1);
full_sos_actual_psd_dim = NaN(num_cases, 1);
full_sos_status = strings(num_cases, 1);
full_sos_log_bound = NaN(num_cases, 1);
full_sos_time = NaN(num_cases, 1);
full_sos_log_gap = NaN(num_cases, 1);
full_sos_flat_order = NaN(num_cases, 1);
full_sos_flat_rank = NaN(num_cases, 1);

row = 0;
for pidx = 1:num_profiles
    weights = config.weight_profiles(pidx, :)';
    for sidx = 1:numel(config.seeds)
        row = row+1;
        profile(row) = config.profile_names(pidx);
        seed(row) = config.seeds(sidx);
        weights_text(row) = sprintf('%.10g,%.10g,%.10g', weights);
        sum_weights(row) = sum(weights);
        integer_weights(row) = all(abs(weights-round(weights)) <= 1e-10);
        fprintf('[%d/%d] profile=%s, seed=%d\n', ...
            row, num_cases, profile(row), seed(row));

        instance = Table_5_6_generate_instance(config.n, config.m, ...
            seed(row), config.epsilon);
        instance.weights = weights;
        local = Table_5_6_solve_sqp(instance, config.num_starts);
        standard = Table_5_6_solve_log_relaxation(instance, false, 1);
        tighter = Table_5_6_solve_log_relaxation(instance, true, config.m+1);

        local_value(row) = local.value;
        local_time(row) = local.time;
        local_successes(row) = local.successful_starts;
        local_solver(row) = local.solver;
        mom_bound(row) = standard.bound;
        mom_time(row) = standard.time;
        mom_log_gap(row) = max(0, standard.bound-local.value);
        mom_gap_per_weight(row) = mom_log_gap(row)/sum(weights);
        mom_multiplicative_gap(row) = exp(min(mom_log_gap(row), log(realmax)));
        lme_bound(row) = tighter.bound;
        lme_time(row) = tighter.time;
        lme_log_gap(row) = max(0, tighter.bound-local.value);
        lme_gap_per_weight(row) = lme_log_gap(row)/sum(weights);
        lme_flat_order(row) = tighter.flat_order;
        lme_flat_rank(row) = tighter.flat_rank;

        if ~integer_weights(row)
            full_sos_status(row) = "not_applicable_noninteger";
            continue
        end
        full_sos_required_order(row) = sum(round(weights));
        full_sos_required_psd_dim(row) = nchoosek(config.n+ ...
            full_sos_required_order(row), full_sos_required_order(row));
        if full_sos_required_psd_dim(row) > config.max_full_sos_psd_dim
            full_sos_status(row) = "skipped_size";
            continue
        end
        try
            full = Table_5_6_solve_direct_product(instance, full_sos_required_order(row), ...
                config.full_sos_max_extra_orders, config.max_full_sos_psd_dim);
            full_sos_status(row) = "solved";
            full_sos_actual_order(row) = full.order;
            full_sos_actual_psd_dim(row) = nchoosek(config.n+full.order, full.order);
            full_sos_log_bound(row) = log(full.bound);
            full_sos_time(row) = full.time;
            full_sos_log_gap(row) = max(0, full_sos_log_bound(row)-local.value);
            full_sos_flat_order(row) = full.flat_order;
            full_sos_flat_rank(row) = full.flat_rank;
        catch err
            full_sos_status(row) = "failed: " + string(err.identifier);
            warning('Full SOS failed for profile %s, seed %d: %s', ...
                profile(row), seed(row), err.message);
        end
    end
end

results = table(profile, seed, n, m, weights_text, sum_weights, ...
    integer_weights, local_value, local_time, local_successes, local_solver, ...
    mom_bound, mom_time, ...
    mom_psd_dim, mom_log_gap, mom_gap_per_weight, mom_multiplicative_gap, ...
    lme_order, lme_psd_dim, lme_bound, lme_time, lme_log_gap, ...
    lme_gap_per_weight, lme_flat_order, lme_flat_rank, ...
    full_sos_required_order, full_sos_required_psd_dim, ...
    full_sos_actual_order, full_sos_actual_psd_dim, full_sos_status, ...
    full_sos_log_bound, full_sos_time, full_sos_log_gap, ...
    full_sos_flat_order, full_sos_flat_rank);
fprintf('\n=== Table 6: MATLAB quantities reported in the paper ===\n');
paper_results = results(:, {'weights_text','mom_log_gap','lme_log_gap', ...
    'full_sos_log_gap','mom_time','lme_time','full_sos_time'});
disp(paper_results)
fprintf('Detailed relaxation results remain available in the variable results.\n');
