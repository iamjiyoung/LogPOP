% scaling results reported in Table 5.


clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

config.m_values = [2 4 6 8 10 12];
config.n_values = [2 4 6 8];
config.fixed_n = 4;
config.fixed_m = 4;
config.seeds = 3;
config.epsilon = 0.05;
config.num_starts = 100;
config.run_lme = true;
config.max_lme_psd_dim = 200;
config.run_full_sos = true;
config.max_full_sos_psd_dim = 250;
config.full_sos_max_extra_orders = 2;

num_m_cases = numel(config.m_values)*numel(config.seeds);
num_n_cases = numel(config.n_values)*numel(config.seeds);
num_cases = num_m_cases+num_n_cases;
study = [repmat("m", num_m_cases, 1); repmat("n", num_n_cases, 1)];
seed = [repelem(config.seeds(:), numel(config.m_values), 1); ...
    repelem(config.seeds(:), numel(config.n_values), 1)];
n = [repmat(config.fixed_n, num_m_cases, 1); ...
    repmat(config.n_values(:), numel(config.seeds), 1)];
m = [repmat(config.m_values(:), numel(config.seeds), 1); ...
    repmat(config.fixed_m, num_n_cases, 1)];
local_value = NaN(num_cases, 1);
local_time = NaN(num_cases, 1);
local_successes = NaN(num_cases, 1);
local_solver = strings(num_cases, 1);
mom_bound = NaN(num_cases, 1);
mom_time = NaN(num_cases, 1);
mom_psd_dim = NaN(num_cases, 1);
mom_full_rank = NaN(num_cases, 1);
mom_flat_order = NaN(num_cases, 1);
mom_flat_rank = NaN(num_cases, 1);
log_gap = NaN(num_cases, 1);
multiplicative_gap = NaN(num_cases, 1);
lme_order = NaN(num_cases, 1);
lme_psd_dim = NaN(num_cases, 1);
lme_status = strings(num_cases, 1);
lme_bound = NaN(num_cases, 1);
lme_time = NaN(num_cases, 1);
lme_log_gap = NaN(num_cases, 1);
lme_full_rank = NaN(num_cases, 1);
lme_flat_order = NaN(num_cases, 1);
lme_flat_rank = NaN(num_cases, 1);
full_sos_order = NaN(num_cases, 1);
full_sos_psd_dim = NaN(num_cases, 1);
full_sos_status = strings(num_cases, 1);
full_sos_product_bound = NaN(num_cases, 1);
full_sos_log_bound = NaN(num_cases, 1);
full_sos_time = NaN(num_cases, 1);
full_sos_log_gap = NaN(num_cases, 1);
full_sos_full_rank = NaN(num_cases, 1);
full_sos_flat_order = NaN(num_cases, 1);
full_sos_flat_rank = NaN(num_cases, 1);

for row = 1:num_cases
    fprintf('[%d/%d] %s-scaling: n=%d, m=%d, seed=%d\n', ...
        row, num_cases, study(row), n(row), m(row), seed(row));

    instance = Table_5_6_generate_instance(n(row), m(row), seed(row), config.epsilon);
    local = Table_5_6_solve_sqp(instance, config.num_starts);
    standard = Table_5_6_solve_log_relaxation(instance, false, 1);

    local_value(row) = local.value;
    local_time(row) = local.time;
    local_successes(row) = local.successful_starts;
    local_solver(row) = local.solver;
    mom_bound(row) = standard.bound;
    mom_time(row) = standard.time;
    mom_psd_dim(row) = nchoosek(n(row)+1, 1);
    mom_full_rank(row) = standard.full_rank;
    mom_flat_order(row) = standard.flat_order;
    mom_flat_rank(row) = standard.flat_rank;
    log_gap(row) = max(0, mom_bound(row)-local_value(row));
    multiplicative_gap(row) = exp(min(log_gap(row), log(realmax)));

    lme_order(row) = m(row)+1;
    lme_psd_dim(row) = nchoosek(n(row)+lme_order(row), lme_order(row));
    if ~config.run_lme
        lme_status(row) = "disabled";
    elseif lme_psd_dim(row) > config.max_lme_psd_dim
        lme_status(row) = "skipped_size";
    else
        try
            tighter = Table_5_6_solve_log_relaxation(instance, true, lme_order(row));
            lme_status(row) = "solved";
            lme_bound(row) = tighter.bound;
            lme_time(row) = tighter.time;
            lme_log_gap(row) = max(0, lme_bound(row)-local_value(row));
            lme_full_rank(row) = tighter.full_rank;
            lme_flat_order(row) = tighter.flat_order;
            lme_flat_rank(row) = tighter.flat_rank;
        catch err
            lme_status(row) = "failed: " + string(err.identifier);
            warning('LME failed for n=%d, m=%d, seed=%d: %s', ...
                n(row), m(row), seed(row), err.message);
        end
    end

    full_sos_order(row) = m(row);
    full_sos_psd_dim(row) = nchoosek(n(row)+full_sos_order(row), full_sos_order(row));
    if ~config.run_full_sos
        full_sos_status(row) = "disabled";
    elseif full_sos_psd_dim(row) > config.max_full_sos_psd_dim
        full_sos_status(row) = "skipped_size";
    else
        try
            full_sos = Table_5_6_solve_direct_product(instance, full_sos_order(row), ...
                config.full_sos_max_extra_orders, config.max_full_sos_psd_dim);
            full_sos_status(row) = "solved";
            full_sos_order(row) = full_sos.order;
            full_sos_psd_dim(row) = nchoosek(n(row)+full_sos.order, full_sos.order);
            full_sos_product_bound(row) = full_sos.bound;
            if full_sos.bound > 0
                full_sos_log_bound(row) = log(full_sos.bound);
                full_sos_log_gap(row) = max(0, ...
                    full_sos_log_bound(row)-local_value(row));
            end
            full_sos_time(row) = full_sos.time;
            full_sos_full_rank(row) = full_sos.full_rank;
            full_sos_flat_order(row) = full_sos.flat_order;
            full_sos_flat_rank(row) = full_sos.flat_rank;
        catch err
            full_sos_status(row) = "failed: " + string(err.identifier);
            warning('Full SOS failed for n=%d, m=%d, seed=%d: %s', ...
                n(row), m(row), seed(row), err.message);
        end
    end
end

results = table(study, seed, n, m, local_value, local_time, local_successes, ...
    local_solver, mom_bound, mom_time, mom_psd_dim, mom_full_rank, ...
    mom_flat_order, mom_flat_rank, log_gap, multiplicative_gap, ...
    lme_order, lme_psd_dim, lme_status, lme_bound, lme_time, lme_log_gap, ...
    lme_full_rank, lme_flat_order, lme_flat_rank, full_sos_order, ...
    full_sos_psd_dim, full_sos_status, full_sos_product_bound, ...
    full_sos_log_bound, full_sos_time, full_sos_log_gap, ...
    full_sos_full_rank, full_sos_flat_order, full_sos_flat_rank);
fprintf('\n=== Table 5: quantities reported in the paper ===\n');
paper_results = results(:, {'study','n','m','log_gap', ...
    'multiplicative_gap','mom_time','lme_psd_dim','full_sos_psd_dim'});
disp(paper_results)
fprintf('Detailed relaxation results remain available in the variable results.\n');
