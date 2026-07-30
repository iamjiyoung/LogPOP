% J. Choi, July 29, 2026
%
function out = Ex_12_iv_fit_lcm_em(patterns, counts, num_starts, seed, options)
%EX_12_IV_FIT_LCM_EM Run the multistart EM calculation for Example 12(iv).
% This helper is called by Ex_12_iv.m and should not be run directly.
% It returns the best log-likelihood, fitted parameters, runtime, and
% diagnostics for all EM starts.
% Fit a two-class binary latent class model by deterministic multistart EM.

if nargin < 3
    num_starts = 100;
end
if nargin < 4
    seed = 1;
end
if nargin < 5
    options = struct;
end
if ~isfield(options, 'max_iterations')
    options.max_iterations = 2000;
end
if ~isfield(options, 'relative_tolerance')
    options.relative_tolerance = 1e-8;
end

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(seed, 'twister');
counts = counts(:);
sample_size = sum(counts);
best_loglik = -Inf;
best_mixing = NaN(1, 2);
best_theta = NaN(2, 4);
final_loglik = NaN(num_starts, 1);
iterations = NaN(num_starts, 1);
started = tic;

for start = 1:num_starts
    mixing = [0.25+0.5*rand, 0];
    mixing(2) = 1-mixing(1);
    theta = 0.15+0.7*rand(2, 4);
    previous_loglik = -Inf;

    for iter = 1:options.max_iterations
        component = component_probabilities(patterns, theta);
        weighted = component.*mixing;
        probabilities = sum(weighted, 2);
        responsibilities = weighted./probabilities;
        weighted_counts = counts.*responsibilities;
        class_counts = sum(weighted_counts, 1);

        mixing = class_counts/sample_size;
        for k = 1:2
            theta(k,:) = (weighted_counts(:,k)'*patterns)/class_counts(k);
        end
        theta = min(max(theta, 1e-10), 1-1e-10);
        mixing = min(max(mixing, 1e-10), 1-1e-10);
        mixing = mixing/sum(mixing);

        component = component_probabilities(patterns, theta);
        probabilities = component*mixing';
        loglik = sum(counts.*log(probabilities));
        if abs(loglik-previous_loglik) <= options.relative_tolerance*(1+abs(loglik))
            break
        end
        previous_loglik = loglik;
    end

    final_loglik(start) = loglik;
    iterations(start) = iter;
    if loglik > best_loglik
        best_loglik = loglik;
        best_mixing = mixing;
        best_theta = theta;
    end
    if start == 1 || mod(start, 10) == 0 || start == num_starts
        fprintf(['EM start %3d/%3d: final log-likelihood = %.10f, ' ...
            'iterations = %d, best = %.10f\n'], ...
            start, num_starts, loglik, iter, best_loglik);
    end
end

rounded_values = round(final_loglik, 6);
out.loglik = best_loglik;
out.mixing = best_mixing;
out.theta = best_theta;
out.time = toc(started);
out.final_loglik = final_loglik;
out.iterations = iterations;
out.num_starts = num_starts;
out.num_distinct_solutions = numel(unique(rounded_values));
out.num_best_solutions = sum(abs(final_loglik-best_loglik) <= 1e-6);
out.options = options;
end

function component = component_probabilities(patterns, theta)
component = ones(size(patterns,1), 2);
for k = 1:2
    for j = 1:4
        component(:,k) = component(:,k).* ...
            theta(k,j).^patterns(:,j).*(1-theta(k,j)).^(1-patterns(:,j));
    end
end
end
