clear
clc

orders = 2:6;
bounds = zeros(size(orders));
flat_orders = NaN(size(orders));
flat_ranks = NaN(size(orders));

for row = 1:length(orders)
    ord = orders(row);

    %%% GloptiPoly code
    mset clear
    mpol x 3

    p{1} = x(1)^2 + 2*x(1)*x(3);
    p{2} = x(2)^2 + 2*x(2)*x(3);
    p{3} = 2*x(1)*x(2);
    p{4} = x(3)^2;

    report_weights = [182 60 17 176];
    solve_weights = report_weights/max(report_weights);

    Ktheta = [1-sum(x) >= 0];

    % A redundant ball constraint makes the description Archimedean.
    Ktheta = [Ktheta; 1-x'*x >= 0];
    for i = 1:length(x)
        Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
    end

    PO = cell(length(p), 1);
    A = cell(length(p), 1);
    bp = cell(length(p), 1);
    c = cell(length(p), 1);
    K = cell(length(p), 1);
    b0 = cell(length(p), 1);
    for i = 1:length(p)
        PO{i} = msdp(max(p{i}), Ktheta, ord);
        [A{i}, bp{i}, c{i}, K{i}, b0{i}] = msedumi(PO{i});
    end

    %%%% YALMIP code
    MCone = K{1};
    mom_y = sdpvar(length(bp{1}), 1);
    c_Aty = c{1}-A{1}'*mom_y;
    MomRelax = [];

    for k = 1:MCone.f
        MomRelax = [MomRelax, c_Aty(k) == 0]; %#ok<AGROW>
    end

    offset = MCone.f;
    for k = 1:MCone.l
        MomRelax = [MomRelax, c_Aty(offset+k) >= 0]; %#ok<AGROW>
    end
    offset = offset+MCone.l;

    moment = [];
    for k = 1:length(MCone.s)
        block_size = MCone.s(k)^2;
        locM = mat(c_Aty(offset+1:offset+block_size));
        if k == 1
            moment = locM;
        end
        MomRelax = [MomRelax, locM >= 0]; %#ok<AGROW>
        offset = offset+block_size;
    end

    obj = 0;
    log_arguments = cell(length(p), 1);
    for i = 1:length(p)
        log_arguments{i} = b0{i}+bp{i}'*mom_y;
        obj = obj-solve_weights(i)*log(log_arguments{i});
    end

    sol = optimize(MomRelax, obj, sdpsettings('solver', 'mosek'));
    if sol.problem ~= 0
        error('MOSEK failed at order %d: %s', ord, sol.info);
    end

    arguments = zeros(1, length(p));
    for i = 1:length(p)
        arguments(i) = value(log_arguments{i});
    end
    bounds(row) = sum(report_weights.*log(arguments));

    % Check rank(M_t) = rank(M_{t-1}) for t = 1,...,k.
    moment_value = value(moment);
    for t = ord:-1:1
        current_size = nchoosek(length(x)+t, t);
        previous_size = nchoosek(length(x)+t-1, t-1);
        current_rank = rank( ...
            moment_value(1:current_size, 1:current_size), 1e-4);
        previous_rank = rank( ...
            moment_value(1:previous_size, 1:previous_size), 1e-4);
        if current_rank == previous_rank
            candidate_moment = ...
                moment_value(1:current_size, 1:current_size);
            candidate_atoms = extractmin(candidate_moment, length(x), t);
            atoms_are_feasible = ~isempty(candidate_atoms);
            for atom_index = 1:length(candidate_atoms)
                atom = candidate_atoms{atom_index}(:);
                atom_violation = max([0; -atom; sum(atom)-1; ...
                    atom'*atom-1]);
                if atom_violation > 1e-4
                    atoms_are_feasible = false;
                    break
                end
            end
            if atoms_are_feasible
                flat_orders(row) = t;
                flat_ranks(row) = current_rank;
                break
            end
        end
    end
end

fprintf('\n=== Standard moment relaxation for the ABO example ===\n');
for row = 1:length(orders)
    fprintf('\nRelaxation order k        : %d\n', orders(row));
    fprintf('Optimal relaxation value  : %.10f\n', bounds(row));
    if isnan(flat_orders(row))
        fprintf('Flat truncation            : not detected\n');
    else
        fprintf('Flat-truncation order t    : %d\n', flat_orders(row));
        fprintf('Rank                       : %d\n', flat_ranks(row));
    end
end
