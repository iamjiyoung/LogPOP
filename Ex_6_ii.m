% five coefficient instances in Example 6(ii).
%
% J. Choi, August 10, 2026

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), ...
    'functions', 'setup_paths.m'))

weights = [
    0.7691 0.6389 0.8931 0.0607 0.1758 0.4163;
    0.1774 0.3959 0.4922 0.4379 0.6354 0.1527;
    0.2920 0.4317 0.0155 0.9841 0.1672 0.1062;
    0.0835 0.6260 0.6609 0.7298 0.8908 0.9823;
    0.3424 0.7360 0.7947 0.5449 0.6862 0.8936];
mom_orders = NaN(5,1);
lme_orders = NaN(5,1);

instance = (1:5)';
rank_mom = NaN(5,1);
flat_t_mom = NaN(5,1);
f_mom = NaN(5,1);
f_lme = NaN(5,1);
rank_lme = NaN(5,1);
flat_t_lme = NaN(5,1);
time_mom = NaN(5,1);
time_lme = NaN(5,1);
f_factorized = NaN(5,1);
time_factorized = NaN(5,1);
maximizer_factorized = NaN(5,3);
atoms_mom = strings(5,1);
atoms_lme = strings(5,1);
atom_values_mom = cell(5,1);
atom_values_lme = cell(5,1);
standard_order_history = cell(5,1);
standard_bound_history = cell(5,1);
standard_rank_history = cell(5,1);
standard_flat_order_history = cell(5,1);
standard_flat_rank_history = cell(5,1);
atom_feasibility_tolerance = 1e-4;

for i = 1:5
    fprintf('\n============================================================\n');
    fprintf('Example 6(ii), instance %d\n', i);
    fprintf('Coefficient vector: %s\n', mat2str(weights(i,:), 10));
    fprintf('============================================================\n');
    standard = [];
    standard_last = [];
    for candidate_order = 2:8
        fprintf('\n--- Standard moment relaxation, order k=%d ---\n', ...
            candidate_order);
        candidate = Ex_6_ii_solve_case(weights(i,:), false, candidate_order);
        if ~isnan(candidate.flat_order)
            max_violation = 0;
            for j = 1:numel(candidate.atoms)
                atom = candidate.atoms{j}(:);
                violation = max([0; -atom; sum(atom)-1; atom'*atom-1]);
                max_violation = max(max_violation, violation);
            end
            if max_violation > atom_feasibility_tolerance
                fprintf(['Rejected numerical flat truncation: ' ...
                    'maximum atom feasibility violation = %.3e.\n'], ...
                    max_violation);
                candidate.flat_order = NaN;
                candidate.flat_rank = NaN;
                candidate.atoms = {};
            end
            if ~isnan(candidate.flat_order) && candidate.flat_rank == 1
                objective_residual = ...
                    abs(candidate.atom_values(1)-candidate.bound);
                if objective_residual > 1e-4
                    fprintf(['Rejected numerical flat truncation: ' ...
                        'objective reconstruction residual = %.3e.\n'], ...
                        objective_residual);
                    candidate.flat_order = NaN;
                    candidate.flat_rank = NaN;
                    candidate.atoms = {};
                end
            end
        end
        standard_last = candidate;
        standard_order_history{i}(end+1) = candidate_order;
        standard_bound_history{i}(end+1) = candidate.bound;
        standard_rank_history{i}(end+1) = candidate.full_rank;
        standard_flat_order_history{i}(end+1) = candidate.flat_order;
        standard_flat_rank_history{i}(end+1) = candidate.flat_rank;
        if ~isnan(candidate.flat_order)
            standard = candidate;
            mom_orders(i) = candidate_order;
            break
        end
    end
    if isempty(standard)
        standard = standard_last;
        fprintf(['No numerically reliable standard flat truncation ' ...
            'was found for instance %d through order 8.\n'], i);
    end
    tighter = [];
    for candidate_order = 3:6
        fprintf('\n--- LME moment relaxation, order k=%d ---\n', ...
            candidate_order);
        candidate = Ex_6_ii_solve_case(weights(i,:), true, candidate_order);
        if ~isnan(candidate.flat_order)
            max_violation = 0;
            for j = 1:numel(candidate.atoms)
                atom = candidate.atoms{j}(:);
                violation = max([0; -atom; sum(atom)-1; atom'*atom-1]);
                max_violation = max(max_violation, violation);
            end
            if max_violation > atom_feasibility_tolerance
                fprintf(['Rejected numerical flat truncation: ' ...
                    'maximum atom feasibility violation = %.3e.\n'], ...
                    max_violation);
                candidate.flat_order = NaN;
                candidate.flat_rank = NaN;
                candidate.atoms = {};
            end
            if ~isnan(candidate.flat_order) && candidate.flat_rank == 1
                objective_residual = ...
                    abs(candidate.atom_values(1)-candidate.bound);
                if objective_residual > 1e-4
                    fprintf(['Rejected numerical flat truncation: ' ...
                        'objective reconstruction residual = %.3e.\n'], ...
                        objective_residual);
                    candidate.flat_order = NaN;
                    candidate.flat_rank = NaN;
                    candidate.atoms = {};
                end
            end
        end
        if ~isnan(candidate.flat_order)
            tighter = candidate;
            lme_orders(i) = candidate_order;
            break
        end
    end
    if isempty(tighter)
        error('No LME flat truncation found for instance %d through order 6.', i);
    end

    fprintf('\n--- Factorized standard moment relaxation, order k=1 ---\n');
    factorized = Ex_6_factorized_solve_case(weights(i,:));
    rank_mom(i) = standard.flat_rank;
    flat_t_mom(i) = standard.flat_order;
    f_mom(i) = standard.bound;
    f_lme(i) = tighter.bound;
    rank_lme(i) = tighter.flat_rank;
    flat_t_lme(i) = tighter.flat_order;
    time_mom(i) = standard.time;
    time_lme(i) = tighter.time;
    f_factorized(i) = factorized.bound;
    time_factorized(i) = factorized.time;
    maximizer_factorized(i,:) = factorized.maximizer;
    atom_text = strings(1, numel(standard.atoms));
    for j = 1:numel(standard.atoms)
        atom_text(j) = string(mat2str(standard.atoms{j}(:)', 10));
    end
    atoms_mom(i) = strjoin(atom_text, '; ');
    atom_values_mom{i} = standard.atom_values;
    atom_text = strings(1, numel(tighter.atoms));
    for j = 1:numel(tighter.atoms)
        atom_text(j) = string(mat2str(tighter.atoms{j}(:)', 10));
    end
    atoms_lme(i) = strjoin(atom_text, '; ');
    atom_values_lme{i} = tighter.atom_values;
end

fprintf('\n=== Example 6(ii): reported results ===\n');
for i = 1:5
    fprintf('\nInstance %d\n', instance(i));
    fprintf('Coefficient vector         : %s\n', mat2str(weights(i,:), 10));
    fprintf('\nStandard moment relaxation\n');
    if isnan(mom_orders(i))
        fprintf('Flat truncation            : not detected through k=8\n');
        fprintf('Reported bound order k     : %d\n', ...
            standard_order_history{i}(end));
    else
        fprintf('Relaxation order k         : %d\n', mom_orders(i));
    end
    fprintf('Optimal relaxation value   : %.10f\n', f_mom(i));
    if ~isnan(flat_t_mom(i))
        fprintf('Flat-truncation order t     : %d\n', flat_t_mom(i));
        fprintf('Rank                        : %d\n', rank_mom(i));
        fprintf('Extracted atoms             : %s\n', atoms_mom(i));
        fprintf('Objectives at atoms         : %s\n', ...
            mat2str(atom_values_mom{i}, 10));
    end
    fprintf('Runtime (seconds)           : %.4f\n', time_mom(i));
    fprintf('\nLME moment relaxation\n');
    fprintf('Relaxation order k         : %d\n', lme_orders(i));
    fprintf('Optimal relaxation value   : %.10f\n', f_lme(i));
    fprintf('Flat-truncation order t     : %d\n', flat_t_lme(i));
    fprintf('Rank                        : %d\n', rank_lme(i));
    fprintf('Extracted atoms             : %s\n', atoms_lme(i));
    fprintf('Objectives at atoms         : %s\n', ...
        mat2str(atom_values_lme{i}, 10));
    fprintf('Runtime (seconds)           : %.4f\n', time_lme(i));
    fprintf('\nFactorized standard moment relaxation\n');
    fprintf('Relaxation order k         : 1\n');
    fprintf('Optimal relaxation value   : %.10f\n', f_factorized(i));
    fprintf('First-moment maximizer     : %s\n', ...
        mat2str(maximizer_factorized(i,:), 10));
    fprintf('Runtime (seconds)           : %.4f\n', time_factorized(i));
end

fprintf('\n=== Example 6(ii), instance 3: standard order history ===\n');
for j = 1:numel(standard_order_history{3})
    if isnan(standard_flat_order_history{3}(j))
        flat_text = 'no';
    else
        flat_text = sprintf('yes (t=%d, rank=%d)', ...
            standard_flat_order_history{3}(j), ...
            standard_flat_rank_history{3}(j));
    end
    fprintf(['k=%d: bound=%.10f, moment rank=%d, ' ...
        'flat truncation=%s\n'], ...
        standard_order_history{3}(j), ...
        standard_bound_history{3}(j), ...
        standard_rank_history{3}(j), flat_text);
end
