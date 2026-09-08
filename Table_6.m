% Compare the log-moment and direct-product Moment-SOS formulations as
% the positive weights change. Run this script directly.

clear
clc

n = 4;
m = 3;
epsilon = 0.05;
seed = 3;
max_product_matrix_size = 250;

weight_profiles = {
    [1; 1; 1]
    [1; 2; 3]
    [2; 5; 10]
    [5; 10; 20]
    [sqrt(2); pi/2; exp(1)/2]
};
weight_labels = [
    "(1,1,1)"
    "(1,2,3)"
    "(2,5,10)"
    "(5,10,20)"
    "(sqrt(2),pi/2,e/2)"
];

% Generate the fixed instance used for the weight comparison.
rng(seed, 'twister');
U = zeros(m, n);
beta = zeros(m, 1);
for j = 1:m
    U(j, :) = randn(1, n);
    U(j, :) = U(j, :)/norm(U(j, :));
    beta(j) = 0.15*(2*rand-1);
end

num_cases = numel(weight_profiles);
log_order = ones(num_cases, 1);
log_matrix_size = (n+1)*ones(num_cases, 1);
log_bound = NaN(num_cases, 1);
log_time = NaN(num_cases, 1);
product_order = NaN(num_cases, 1);
product_matrix_size = NaN(num_cases, 1);
product_log_bound = NaN(num_cases, 1);
product_time = NaN(num_cases, 1);
product_status = strings(num_cases, 1);

for row = 1:num_cases
    weights = weight_profiles{row};
    fprintf('\n[%d/%d] weights = %s\n', row, num_cases, weight_labels(row));

    mset clear
    mpol x 4

    p = cell(m, 1);
    for j = 1:m
        affine = beta(j);
        for i = 1:n
            affine = affine+U(j, i)*x(i);
        end
        p{j} = epsilon+affine^2;
    end

    Ktheta = [x'*x <= 1; 1-sum(x) >= 0];
    for i = 1:n
        Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
    end

    % Standard order-one log-moment relaxation.
    for j = 1:m
        problem{j} = msdp(max(p{j}), Ktheta, 1); %#ok<SAGROW>
        [A{j}, bp{j}, cc{j}, cone{j}, b0{j}] = msedumi(problem{j}); %#ok<SAGROW>
    end

    mom_y = sdpvar(length(bp{1}), 1);
    c_Aty = cc{1}-A{1}'*mom_y;
    constraints = [];
    for i = 1:cone{1}.f
        constraints = [constraints, c_Aty(i) == 0]; %#ok<AGROW>
    end
    offset = cone{1}.f;
    for i = 1:cone{1}.l
        constraints = [constraints, c_Aty(offset+i) >= 0]; %#ok<AGROW>
    end
    offset = offset+cone{1}.l;
    for i = 1:numel(cone{1}.s)
        block_size = cone{1}.s(i)^2;
        locM = mat(c_Aty(offset+1:offset+block_size));
        constraints = [constraints, locM >= 0]; %#ok<AGROW>
        offset = offset+block_size;
    end

    objective = 0;
    for j = 1:m
        objective = objective-weights(j)*log(b0{j}+bp{j}'*mom_y);
    end
    settings = sdpsettings('solver', 'mosek');
    optimize(constraints, objective, settings); % Unreported warm-up solve
    started = tic;
    solution = optimize(constraints, objective, settings);
    log_time(row) = toc(started);
    if solution.problem ~= 0
        error('The log-moment relaxation failed: %s', solution.info);
    end
    log_bound(row) = -value(objective);

    % The direct product is a polynomial only for positive integer weights.
    integer_weights = round(weights);
    if any(abs(weights-integer_weights) > 1e-10)
        product_status(row) = "not applicable";
        continue
    end

    product_order(row) = sum(integer_weights);
    product_matrix_size(row) = nchoosek(n+product_order(row), product_order(row));
    if product_matrix_size(row) > max_product_matrix_size
        product_status(row) = "not launched";
        continue
    end

    product_objective = 1;
    for j = 1:m
        product_objective = product_objective*p{j}^integer_weights(j);
    end
    product_problem = msdp(max(product_objective), Ktheta, product_order(row));
    [product_A, product_bp, product_c, product_cone, product_b0] = ...
        msedumi(product_problem);

    product_y = sdpvar(length(product_bp), 1);
    product_c_Aty = product_c-product_A'*product_y;
    product_constraints = [];
    for i = 1:product_cone.f
        product_constraints = [product_constraints, ...
            product_c_Aty(i) == 0]; %#ok<AGROW>
    end
    offset = product_cone.f;
    for i = 1:product_cone.l
        product_constraints = [product_constraints, ...
            product_c_Aty(offset+i) >= 0]; %#ok<AGROW>
    end
    offset = offset+product_cone.l;
    for i = 1:numel(product_cone.s)
        block_size = product_cone.s(i)^2;
        locM = mat(product_c_Aty(offset+1:offset+block_size));
        product_constraints = [product_constraints, locM >= 0]; %#ok<AGROW>
        offset = offset+block_size;
    end

    linear_objective = -(product_b0+product_bp'*product_y);
    optimize(product_constraints, linear_objective, settings); % Warm-up
    started = tic;
    solution = optimize(product_constraints, linear_objective, settings);
    product_time(row) = toc(started);
    if solution.problem ~= 0
        product_status(row) = "failed";
    else
        product_status(row) = "solved";
        product_bound = value(product_b0+product_bp'*product_y);
        if product_bound > 0
            product_log_bound(row) = log(product_bound);
        end
    end
end

results = table(weight_labels, log_order, log_matrix_size, log_bound, ...
    log_time, product_order, product_matrix_size, product_log_bound, ...
    product_time, product_status);

fprintf('\n=== Weight experiment ===\n');
disp(results)
