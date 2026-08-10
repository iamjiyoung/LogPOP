% EM comparisons for Examples 4, 11, and 12.
%
% The script prints the best objective value, returned parameter vector,
% number of iterations, and wall-clock runtime.  EM returns feasible
% candidates but does not by itself provide a global-optimality certificate.
%
% J. Choi, August 5, 2026

clear
clc

max_iterations = 5000;
relative_tolerance = 1e-10;

%% Example 4: ABO phenotype likelihood
counts = [182, 60, 17, 176];
sample_size = sum(counts);
num_starts = 100;
rng(20411004, 'twister');
best_value = -Inf;
best_x = NaN(1, 3);
best_iterations = NaN;
started = tic;

for start = 1:num_starts
    x = rand(1, 3);
    x = x/sum(x);
    previous_value = -Inf;
    for iteration = 1:max_iterations
        pA = x(1)^2+2*x(1)*x(3);
        pB = x(2)^2+2*x(2)*x(3);
        expected_AA = counts(1)*x(1)^2/pA;
        expected_AO = counts(1)*2*x(1)*x(3)/pA;
        expected_BB = counts(2)*x(2)^2/pB;
        expected_BO = counts(2)*2*x(2)*x(3)/pB;
        x_new = [
            2*expected_AA+expected_AO+counts(3), ...
            2*expected_BB+expected_BO+counts(3), ...
            expected_AO+expected_BO+2*counts(4)]/(2*sample_size);
        x = x_new/sum(x_new);
        probabilities = [
            x(1)^2+2*x(1)*x(3), ...
            x(2)^2+2*x(2)*x(3), ...
            2*x(1)*x(2), ...
            x(3)^2];
        value = sum(counts.*log(probabilities));
        if abs(value-previous_value) <= relative_tolerance*(1+abs(value))
            break
        end
        previous_value = value;
    end
    if value > best_value
        best_value = value;
        best_x = x;
        best_iterations = iteration;
    end
end
runtime = toc(started);

fprintf('\n=== Example 4: 100-start ABO EM ===\n');
fprintf('Best feasible objective       : %.10f\n', best_value);
fprintf('Returned optimizer candidate  : %s\n', mat2str(best_x, 10));
fprintf('Iterations for best start     : %d\n', best_iterations);
fprintf('Total EM runtime (seconds)    : %.6f\n', runtime);
fprintf('Globality certified by EM     : no\n');

%% Example 11: paternity mixture weights
Ndata = {
    [77 23], [63 37], [49 40 11], [83 2 15], [63 17 20], ...
    [59 8 16 17], [7 9 4 33 47], [39 38 23], [29 21 88 62]};
Pdata = {
    [0.5 1; 0.5 0], ...
    [0.5 1; 0.5 0], ...
    [0.5 0.875; 0.25 0.125; 0.25 0], ...
    [0.5 0.25; 0.5 0.5; 0 0.25], ...
    [0.5 0.25; 0.5 0.5; 0 0.25], ...
    [0.25 0.5; 0.25 0.5; 0.25 0; 0.25 0], ...
    [0.25 0; 0.25 0; 0.25 0; 0.25 0.5; 0 0.5], ...
    [0.5 0 0.875; 0.25 0.75 0.125; 0.25 0.25 0], ...
    [0.5 0 0; 0.25 0.75 0.25; 0.25 0.25 0.5; 0 0 0.25]};

fprintf('\n=== Example 11: paternity EM ===\n');
for instance = 1:numel(Ndata)
    counts = Ndata{instance}(:);
    P = Pdata{instance};
    sample_size = sum(counts);
    x = ones(1, size(P, 2))/size(P, 2);
    previous_value = -Inf;
    started = tic;
    for iteration = 1:max_iterations
        probabilities = P*x';
        responsibilities = P.*x./probabilities;
        expected_counts = counts.*responsibilities;
        x = sum(expected_counts, 1)/sample_size;
        probabilities = P*x';
        value = sum(counts.*log(probabilities));
        if abs(value-previous_value) <= relative_tolerance*(1+abs(value))
            break
        end
        previous_value = value;
    end
    runtime = toc(started);
    fprintf('\nInstance %d\n', instance);
    fprintf('Best feasible objective       : %.10f\n', value);
    fprintf('Returned optimizer candidate  : %s\n', mat2str(x, 10));
    fprintf('Iterations                    : %d\n', iteration);
    fprintf('EM runtime (seconds)          : %.6f\n', runtime);
    fprintf('Globality certified by EM     : no\n');
end

%% Example 12: latent class models
case_names = ["i", "ii", "iii", "iv"];
d_values = [1, 1, 2, 4];
class_values = [2, 3, 2, 2];
count_data = {
    [303 197], ...
    [227 273], ...
    [131 174 49 146], ...
    [30 26 39 26 27 19 28 30 24 25 39 33 45 31 35 43]};
num_starts = 100;

fprintf('\n=== Example 12: 100-start LCM EM ===\n');
for case_index = 1:4
    d = d_values(case_index);
    number_of_classes = class_values(case_index);
    counts = count_data{case_index}(:);
    patterns = dec2bin(0:(2^d-1))-'0';
    sample_size = sum(counts);
    rng(20361000+case_index, 'twister');
    best_value = -Inf;
    best_mixing = NaN(1, number_of_classes);
    best_theta = NaN(number_of_classes, d);
    best_iterations = NaN;
    started = tic;

    for start = 1:num_starts
        if number_of_classes == 2
            % Match the initialization used by Ex_12_iv.m.
            mixing = [0.25+0.5*rand, 0];
            mixing(2) = 1-mixing(1);
        else
            mixing = rand(1, number_of_classes);
            mixing = mixing/sum(mixing);
        end
        theta = 0.15+0.7*rand(number_of_classes, d);
        previous_value = -Inf;
        for iteration = 1:max_iterations
            component = ones(size(patterns, 1), number_of_classes);
            for class = 1:number_of_classes
                for variable = 1:d
                    component(:, class) = component(:, class).* ...
                        theta(class, variable).^patterns(:, variable).* ...
                        (1-theta(class, variable)).^(1-patterns(:, variable));
                end
            end
            weighted = component.*mixing;
            probabilities = sum(weighted, 2);
            responsibilities = weighted./probabilities;
            weighted_counts = counts.*responsibilities;
            class_counts = sum(weighted_counts, 1);
            mixing = class_counts/sample_size;
            for class = 1:number_of_classes
                theta(class, :) = ...
                    (weighted_counts(:, class)'*patterns)/class_counts(class);
            end
            theta = min(max(theta, 1e-10), 1-1e-10);
            mixing = min(max(mixing, 1e-10), 1-1e-10);
            mixing = mixing/sum(mixing);

            component = ones(size(patterns, 1), number_of_classes);
            for class = 1:number_of_classes
                for variable = 1:d
                    component(:, class) = component(:, class).* ...
                        theta(class, variable).^patterns(:, variable).* ...
                        (1-theta(class, variable)).^(1-patterns(:, variable));
                end
            end
            probabilities = component*mixing';
            value = sum(counts.*log(probabilities));
            if abs(value-previous_value) <= relative_tolerance*(1+abs(value))
                break
            end
            previous_value = value;
        end
        if value > best_value
            best_value = value;
            best_mixing = mixing;
            best_theta = theta;
            best_iterations = iteration;
        end
    end
    runtime = toc(started);

    fprintf('\nCase (%s)\n', case_names(case_index));
    fprintf('Best feasible objective       : %.10f\n', best_value);
    fprintf('Returned mixing weights       : %s\n', mat2str(best_mixing, 10));
    fprintf('Returned class probabilities  : %s\n', mat2str(best_theta, 10));
    fprintf('Iterations for best start     : %d\n', best_iterations);
    fprintf('Total EM runtime (seconds)    : %.6f\n', runtime);
    fprintf('Globality certified by EM     : no\n');
end
