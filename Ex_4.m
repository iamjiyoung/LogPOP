% Example 4: ABO blood-group likelihood.

clear
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))

orders = (2:6)';
mom_bound = zeros(size(orders));
lme_bound = zeros(size(orders));
mom_time = zeros(size(orders));
lme_time = zeros(size(orders));
mom_flat_t = NaN(size(orders));
mom_flat_r = NaN(size(orders));
lme_flat_t = NaN(size(orders));
lme_flat_r = NaN(size(orders));
lme_atoms = cell(size(orders));

for row = 1:numel(orders)
    ord = orders(row);

    % method = 1: standard moment relaxation
    % method = 2: LME moment relaxation
    for method = 1:2
        use_lme = method == 2;

        %% Define the polynomial optimization problem
        mset clear
        mpol('x', 3)

        p{1} = x(1)^2+2*x(1)*x(3);
        p{2} = x(2)^2+2*x(2)*x(3);
        p{3} = 2*x(1)*x(2);
        p{4} = x(3)^2;

        report_weights = [182 60 17 176];
        solve_weights = report_weights/max(report_weights);

        Ktheta = [1-sum(x) >= 0];
        for i = 1:length(x)
            Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
        end

        %% Add the LME equations when requested
        flat_degree = 1;
        if use_lme
            flat_degree = 2;
            lme0 = 2*sum(solve_weights);
            lme(1) = lme0*x(1)*(x(1)+2*x(3)) ...
                -(solve_weights(1)+solve_weights(3))*(x(1)+2*x(3)) ...
                -solve_weights(1)*x(1);
            lme(2) = lme0*x(2)*(x(2)+2*x(3)) ...
                -(solve_weights(2)+solve_weights(3))*(x(2)+2*x(3)) ...
                -solve_weights(2)*x(2);
            lme(3) = lme0*x(3)*(x(1)+2*x(3))*(x(2)+2*x(3)) ...
                -2*solve_weights(4)*(x(1)+2*x(3))*(x(2)+2*x(3)) ...
                -2*solve_weights(1)*x(3)*(x(2)+2*x(3)) ...
                -2*solve_weights(2)*x(3)*(x(1)+2*x(3));
            for i = 1:length(x)
                Ktheta = [Ktheta; lme(i) >= 0; ...
                    lme(i)*x(i) == 0]; %#ok<AGROW>
            end
            Ktheta = [Ktheta; lme0*(1-sum(x)) == 0];
        end

        %% Build the moment relaxation with GloptiPoly
        num_terms = numel(p);
        PO = cell(num_terms, 1);
        A = cell(num_terms, 1);
        bp = cell(num_terms, 1);
        c = cell(num_terms, 1);
        K = cell(num_terms, 1);
        b0 = cell(num_terms, 1);
        for i = 1:num_terms
            PO{i} = msdp(max(p{i}), Ktheta, ord);
            [A{i}, bp{i}, c{i}, K{i}, b0{i}] = msedumi(PO{i});
        end

        %% Write the equality, inequality, and PSD moment constraints
        moment_cone = K{1};
        y = sdpvar(length(bp{1}), 1);
        cone_coordinates = c{1}-A{1}'*y;
        constraints = [];

        for i = 1:moment_cone.f
            constraints = [constraints, ...
                cone_coordinates(i) == 0]; %#ok<AGROW>
        end

        offset = moment_cone.f;
        for i = 1:moment_cone.l
            constraints = [constraints, ...
                cone_coordinates(offset+i) >= 0]; %#ok<AGROW>
        end
        offset = offset+moment_cone.l;

        moment_matrix = [];
        for i = 1:numel(moment_cone.s)
            block_size = moment_cone.s(i)^2;
            psd_block = mat(cone_coordinates(offset+1:offset+block_size));
            constraints = [constraints, psd_block >= 0]; %#ok<AGROW>
            if i == 1
                moment_matrix = psd_block;
            end
            offset = offset+block_size;
        end

        %% Form and solve the logarithmic objective
        log_arguments = cell(num_terms, 1);
        objective = 0;
        for i = 1:num_terms
            log_arguments{i} = b0{i}+bp{i}'*y;
            objective = objective-solve_weights(i)*log(log_arguments{i});
        end

        started = tic;
        solution = optimize(constraints, objective, ...
            sdpsettings('solver', 'mosek'));
        elapsed = toc(started);
        if solution.problem ~= 0
            error('MOSEK failed: %s', solution.info);
        end

        arguments = zeros(1, num_terms);
        for i = 1:num_terms
            arguments(i) = value(log_arguments{i});
        end
        bound = sum(report_weights.*log(arguments));
        moment = value(moment_matrix);

        %% Check rank(M_t) = rank(M_{t-flat_degree})
        flat_t = NaN;
        flat_r = NaN;
        atoms = {};
        previous_rng = rng;
        rng(0, 'twister');
        for candidate = ord:-1:flat_degree
            current_size = nchoosek(length(x)+candidate, candidate);
            previous_order = candidate-flat_degree;
            previous_size = nchoosek(length(x)+previous_order, previous_order);
            current_moment = moment(1:current_size, 1:current_size);
            current_rank = rank(current_moment, 1e-4);
            previous_rank = rank( ...
                moment(1:previous_size, 1:previous_size), 1e-4);
            if current_rank == previous_rank
                flat_t = candidate;
                flat_r = current_rank;
                atoms = extractmin(current_moment, length(x), candidate);
                break
            end
        end
        rng(previous_rng);

        %% Store the result
        if use_lme
            lme_bound(row) = bound;
            lme_time(row) = elapsed;
            lme_flat_t(row) = flat_t;
            lme_flat_r(row) = flat_r;
            lme_atoms{row} = atoms;
        else
            mom_bound(row) = bound;
            mom_time(row) = elapsed;
            mom_flat_t(row) = flat_t;
            mom_flat_r(row) = flat_r;
        end
    end
end

results = table(orders, mom_bound, lme_bound, mom_time, lme_time, ...
    mom_flat_t, mom_flat_r, lme_flat_t, lme_flat_r, 'VariableNames', ...
    {'k','f_mom','f_lme','time_mom','time_lme', ...
    'flat_t_mom','rank_mom','flat_t_lme','rank_lme'});
fprintf('\n=== Example 4: reported results ===\n');

for row = 1:numel(orders)
    fprintf('\nStandard moment relaxation\n');
    fprintf('Relaxation order k        : %d\n', orders(row));
    fprintf('Optimal relaxation value  : %.10f\n', mom_bound(row));
    if isnan(mom_flat_t(row))
        fprintf('Flat truncation            : not detected\n');
    else
        fprintf('Flat-truncation order t    : %d\n', mom_flat_t(row));
        fprintf('Rank                       : %d\n', mom_flat_r(row));
    end
    fprintf('Runtime (seconds)          : %.4f\n', mom_time(row));

    fprintf('\nLME moment relaxation\n');
    fprintf('Relaxation order k        : %d\n', orders(row));
    fprintf('Optimal relaxation value  : %.10f\n', lme_bound(row));
    if isnan(lme_flat_t(row))
        fprintf('Flat truncation            : not detected\n');
    else
        fprintf('Flat-truncation order t    : %d\n', lme_flat_t(row));
        fprintf('Rank                       : %d\n', lme_flat_r(row));
        for atom_index = 1:numel(lme_atoms{row})
            fprintf('Extracted optimizer %d      : %s\n', atom_index, ...
                mat2str(lme_atoms{row}{atom_index}(:)', 10));
        end
    end
    fprintf('Runtime (seconds)          : %.4f\n', lme_time(row));
end
