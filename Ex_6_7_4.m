% Non-saturated latent class model in Example 6.7(iv).
%
% This script solves the order-three moment relaxation and performs a
% 100-start EM calculation for the fixed data set reported in the paper.

clear
clc
mset clear
mset('verbose', true)

%% Fixed data and the saturated likelihood

patterns = dec2bin(0:15)-'0';
counts = [30; 26; 39; 26; 27; 19; 28; 30; ...
          24; 25; 39; 33; 45; 31; 35; 43];
sample_size = sum(counts);

C_obs = [30 26 39 26; ...
         27 19 28 30; ...
         24 25 39 33; ...
         45 31 35 43];
saturated_loglik = sum(counts.*log(counts/sample_size));

%% Order-three moment relaxation

% x(1) is the first mixing weight. The second mixing weight is 1-x(1).
% x(2:5) and x(6:9) are the four Bernoulli probabilities in the two
% latent classes, respectively.
mpol x 9

p = cell(16,1);
for outcome = 1:16
    component_1 = 1;
    component_2 = 1;
    for variable = 1:4
        component_1 = component_1*x(1+variable)^patterns(outcome,variable)* ...
            (1-x(1+variable))^(1-patterns(outcome,variable));
        component_2 = component_2*x(5+variable)^patterns(outcome,variable)* ...
            (1-x(5+variable))^(1-patterns(outcome,variable));
    end
    p{outcome} = x(1)*component_1+(1-x(1))*component_2;
end

Ktheta = 9-x'*x >= 0;
for variable = 1:9
    Ktheta = [Ktheta; x(variable) >= 0; 1-x(variable) >= 0]; %#ok<AGROW>
end

order = 3;
solve_weights = counts/max(counts);
PO = cell(16,1);
A = cell(16,1);
bp = cell(16,1);
c = cell(16,1);
K = cell(16,1);
b0 = cell(16,1);

for outcome = 1:16
    PO{outcome} = msdp(max(p{outcome}),Ktheta,order);
    [A{outcome},bp{outcome},c{outcome},K{outcome},b0{outcome}] = ...
        msedumi(PO{outcome});
end

moment_cone = K{1};
y = sdpvar(length(bp{1}),1);
cone_coordinates = c{1}-A{1}'*y;
constraints = [];

for index = 1:moment_cone.f
    constraints = [constraints,cone_coordinates(index) == 0]; %#ok<AGROW>
end

offset = moment_cone.f;
for index = 1:moment_cone.l
    constraints = [constraints,cone_coordinates(offset+index) >= 0]; %#ok<AGROW>
end
offset = offset+moment_cone.l;

for block = 1:length(moment_cone.s)
    block_size = moment_cone.s(block)^2;
    localizing_matrix = mat(cone_coordinates(offset+1:offset+block_size));
    constraints = [constraints,localizing_matrix >= 0]; %#ok<AGROW>
    offset = offset+block_size;
end

log_arguments = cell(16,1);
objective = 0;
for outcome = 1:16
    log_arguments{outcome} = b0{outcome}+bp{outcome}'*y;
    objective = objective-solve_weights(outcome)*log(log_arguments{outcome});
end

settings = sdpsettings('solver','mosek','verbose',1);
warmup = optimize(constraints,objective,settings);
if warmup.problem ~= 0
    error('MOSEK warm-up failed: %s',warmup.info)
end

started = tic;
solution = optimize(constraints,objective,settings);
moment_time = toc(started);
if solution.problem ~= 0
    error('MOSEK failed: %s',solution.info)
end

arguments = zeros(16,1);
for outcome = 1:16
    arguments(outcome) = value(log_arguments{outcome});
end
moment_bound = sum(counts.*log(arguments));

%% A 100-start EM calculation

number_of_starts = 100;
maximum_iterations = 5000;
relative_tolerance = 1e-10;
rng(20361004,'twister')

best_em_loglik = -Inf;
best_mixing = NaN(1,2);
best_theta = NaN(2,4);
final_loglik = NaN(number_of_starts,1);
iterations = NaN(number_of_starts,1);

started = tic;
for start = 1:number_of_starts
    mixing = [0.25+0.5*rand,0];
    mixing(2) = 1-mixing(1);
    theta = 0.15+0.7*rand(2,4);
    previous_loglik = -Inf;

    for iteration = 1:maximum_iterations
        component = ones(16,2);
        for latent_class = 1:2
            for variable = 1:4
                component(:,latent_class) = component(:,latent_class).* ...
                    theta(latent_class,variable).^patterns(:,variable).* ...
                    (1-theta(latent_class,variable)).^(1-patterns(:,variable));
            end
        end

        weighted = component.*mixing;
        probabilities = sum(weighted,2);
        responsibilities = weighted./probabilities;
        weighted_counts = counts.*responsibilities;
        class_counts = sum(weighted_counts,1);

        mixing = class_counts/sample_size;
        for latent_class = 1:2
            theta(latent_class,:) = ...
                (weighted_counts(:,latent_class)'*patterns)/class_counts(latent_class);
        end
        theta = min(max(theta,1e-10),1-1e-10);
        mixing = min(max(mixing,1e-10),1-1e-10);
        mixing = mixing/sum(mixing);

        component = ones(16,2);
        for latent_class = 1:2
            for variable = 1:4
                component(:,latent_class) = component(:,latent_class).* ...
                    theta(latent_class,variable).^patterns(:,variable).* ...
                    (1-theta(latent_class,variable)).^(1-patterns(:,variable));
            end
        end
        probabilities = component*mixing';
        loglik = sum(counts.*log(probabilities));

        if abs(loglik-previous_loglik) <= ...
                relative_tolerance*(1+abs(loglik))
            break
        end
        previous_loglik = loglik;
    end

    final_loglik(start) = loglik;
    iterations(start) = iteration;
    if loglik > best_em_loglik
        best_em_loglik = loglik;
        best_mixing = mixing;
        best_theta = theta;
    end

    if start == 1 || mod(start,10) == 0
        fprintf(['EM start %3d/%3d: log-likelihood = %.10f, ' ...
            'iterations = %d, best = %.10f\n'], ...
            start,number_of_starts,loglik,iteration,best_em_loglik)
    end
end
em_time = toc(started);

%% Quantities reported in the paper

fprintf('\n=== Example 6.7(iv): reported results ===\n')
fprintf('Observed variables d             : 4\n')
fprintf('Latent classes T                 : 2\n')
fprintf('Sample size                      : %d\n',sample_size)
fprintf('Outcome counts                   : %s\n',mat2str(counts'))
fprintf('Determinant of count flattening  : %.0f\n',det(C_obs))
fprintf('Saturated log-likelihood         : %.10f\n',saturated_loglik)
fprintf('Order-three moment upper bound   : %.10f\n',moment_bound)
fprintf('Best 100-start EM value          : %.10f\n',best_em_loglik)
fprintf('Moment bound minus EM value      : %.10f\n',moment_bound-best_em_loglik)
fprintf('Moment runtime (seconds)         : %.4f\n',moment_time)
fprintf('100-start EM runtime (seconds)   : %.4f\n',em_time)
fprintf('Best mixing weights              : %s\n',mat2str(best_mixing,8))
fprintf('Best Bernoulli probabilities     :\n')
disp(best_theta)
fprintf('EM starts reaching iteration cap : %d\n', ...
    sum(iterations == maximum_iterations))
