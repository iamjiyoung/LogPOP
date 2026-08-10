function out = Table_5_6_solve_sqp(instance, num_starts)
%TABLE_6_7_SOLVE_SQP Compute the best feasible value by multistart SQP.
% This helper is called by Table_5.m and Table_6.m and should not be run
% directly. It solves the original nonlinear problem from deterministic
% starts and returns the best feasible value, point, runtime, and status.
% Find a reproducible feasible lower bound by multistart local optimization.
%
% J. Choi, July 29, 2026

if nargin < 2
    num_starts = 20;
end

n = instance.n;
starts = zeros(n, 0);
starts(:, end+1) = zeros(n, 1);
starts(:, end+1) = ones(n, 1)/(n+1);
starts = [starts, eye(n)]; %#ok<AGROW>

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(100000 + instance.seed, 'twister');
while size(starts, 2) < num_starts
    y = -log(rand(n+1, 1));
    y = y/sum(y);
    starts(:, end+1) = y(1:n); %#ok<AGROW>
end
starts = starts(:, 1:num_starts);

objective = @(x) -sum(instance.weights .* ...
    log(instance.epsilon + (instance.U*x + instance.beta).^2));
use_fmincon = exist('fmincon', 'file') == 2;
if use_fmincon
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'OptimalityTolerance', 1e-9, 'StepTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10, 'MaxIterations', 1000);
    solver_name = 'fmincon-sqp';
else
    options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-9, ...
        'MaxIter', 5000, 'MaxFunEvals', 20000);
    solver_name = 'fminsearch-softmax';
end

best_value = -Inf;
best_x = NaN(n, 1);
successful_starts = 0;
started = tic;
for i = 1:size(starts, 2)
    start_value = -objective(starts(:, i));
    if start_value > best_value
        best_value = start_value;
        best_x = starts(:, i);
    end
    if use_fmincon
        [candidate_x, negative_value, exitflag] = fmincon(objective, starts(:, i), ...
            ones(1, n), 1, [], [], zeros(n, 1), [], [], options);
    else
        initial_y = simplex_to_softmax(starts(:, i));
        transformed_objective = @(y) objective(softmax_to_simplex(y));
        [candidate_y, negative_value, exitflag] = fminsearch( ...
            transformed_objective, initial_y, options);
        candidate_x = softmax_to_simplex(candidate_y);
    end
    if exitflag > 0
        successful_starts = successful_starts + 1;
    end
    candidate_value = -negative_value;
    if all(candidate_x >= -1e-7) && sum(candidate_x) <= 1+1e-7 ...
            && candidate_value > best_value
        best_value = candidate_value;
        best_x = candidate_x;
    end
end

out.value = best_value;
out.x = best_x;
out.time = toc(started);
out.num_starts = size(starts, 2);
out.successful_starts = successful_starts;
out.solver = solver_name;
end

function x = softmax_to_simplex(y)
% The omitted softmax coordinate is the simplex slack 1-sum(x).
shift = max([0; y(:)]);
weights = exp([y(:); 0]-shift);
weights = weights/sum(weights);
x = weights(1:end-1);
end

function y = simplex_to_softmax(x)
floor_value = 1e-8;
x = max(x(:), floor_value);
slack = max(1-sum(x), floor_value);
y = log(x/slack);
end
