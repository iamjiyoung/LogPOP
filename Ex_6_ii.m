% five coefficient instances in Example 6(ii).
%
% J. Choi, July 31, 2026

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
atoms_mom = strings(5,1);
atoms_lme = strings(5,1);

for i = 1:5
    fprintf('\n============================================================\n');
    fprintf('Example 6(ii), instance %d\n', i);
    fprintf('Coefficient vector: %s\n', mat2str(weights(i,:), 10));
    fprintf('============================================================\n');
    standard = [];
    for candidate_order = 2:8
        fprintf('\n--- Standard moment relaxation, order k=%d ---\n', ...
            candidate_order);
        candidate = Ex_6_ii_solve_case(weights(i,:), false, candidate_order);
        if ~isnan(candidate.flat_order)
            standard = candidate;
            mom_orders(i) = candidate_order;
            break
        end
    end
    if isempty(standard)
        error('No flat truncation found for instance %d through order 8.', i);
    end
    tighter = [];
    for candidate_order = 3:6
        fprintf('\n--- LME moment relaxation, order k=%d ---\n', ...
            candidate_order);
        candidate = Ex_6_ii_solve_case(weights(i,:), true, candidate_order);
        if ~isnan(candidate.flat_order)
            tighter = candidate;
            lme_orders(i) = candidate_order;
            break
        end
    end
    if isempty(tighter)
        error('No LME flat truncation found for instance %d through order 6.', i);
    end
    rank_mom(i) = standard.flat_rank;
    flat_t_mom(i) = standard.flat_order;
    f_mom(i) = standard.bound;
    f_lme(i) = tighter.bound;
    rank_lme(i) = tighter.flat_rank;
    flat_t_lme(i) = tighter.flat_order;
    time_mom(i) = standard.time;
    time_lme(i) = tighter.time;
    atom_text = strings(1, numel(standard.atoms));
    for j = 1:numel(standard.atoms)
        atom_text(j) = string(mat2str(standard.atoms{j}(:)', 10));
    end
    atoms_mom(i) = strjoin(atom_text, '; ');
    atom_text = strings(1, numel(tighter.atoms));
    for j = 1:numel(tighter.atoms)
        atom_text(j) = string(mat2str(tighter.atoms{j}(:)', 10));
    end
    atoms_lme(i) = strjoin(atom_text, '; ');
end

fprintf('\n=== Example 6(ii): reported results ===\n');
for i = 1:5
    fprintf('\nInstance %d\n', instance(i));
    fprintf('Coefficient vector         : %s\n', mat2str(weights(i,:), 10));
    fprintf('\nStandard moment relaxation\n');
    fprintf('Relaxation order k         : %d\n', mom_orders(i));
    fprintf('Optimal relaxation value   : %.10f\n', f_mom(i));
    fprintf('Flat-truncation order t     : %d\n', flat_t_mom(i));
    fprintf('Rank                        : %d\n', rank_mom(i));
    fprintf('Extracted atoms             : %s\n', atoms_mom(i));
    fprintf('Runtime (seconds)           : %.4f\n', time_mom(i));
    fprintf('\nLME moment relaxation\n');
    fprintf('Relaxation order k         : %d\n', lme_orders(i));
    fprintf('Optimal relaxation value   : %.10f\n', f_lme(i));
    fprintf('Flat-truncation order t     : %d\n', flat_t_lme(i));
    fprintf('Rank                        : %d\n', rank_lme(i));
    fprintf('Extracted atoms             : %s\n', atoms_lme(i));
    fprintf('Runtime (seconds)           : %.4f\n', time_lme(i));
end
