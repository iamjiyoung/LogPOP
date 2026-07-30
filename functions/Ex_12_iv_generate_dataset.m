function data = Ex_12_iv_generate_dataset(sample_size, seed)
%EX_12_IV_GENERATE_DATASET Generate the fixed non-saturated LCM dataset.
% This helper is called by Ex_12_iv.m and should not be run directly.
% It returns the binary response patterns, fixed counts, generating
% parameters, and associated cell probabilities.
% Generate a fixed weak-separation binary LCM data set with d=4 and K=2.

if nargin < 1
    sample_size = 2000;
end
if nargin < 2
    seed = 20260710;
end

patterns = dec2bin(0:15)-'0';
mixing = [0.5 0.5];
theta = [0.38 0.42 0.46 0.50; 0.62 0.58 0.54 0.50];
probabilities = lcm_probabilities(patterns, mixing, theta);

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(seed, 'twister');
draws = rand(sample_size, 1);
edges = [0; cumsum(probabilities(:))];
counts = zeros(16, 1);
for r = 1:16
    counts(r) = sum(draws >= edges(r) & draws < edges(r+1));
end

data.patterns = patterns;
data.counts = counts;
data.sample_size = sample_size;
data.seed = seed;
data.true_mixing = mixing;
data.true_theta = theta;
data.true_probabilities = probabilities;
end

function probabilities = lcm_probabilities(patterns, mixing, theta)
num_patterns = size(patterns, 1);
probabilities = zeros(num_patterns, 1);
for r = 1:num_patterns
    for k = 1:2
        component = 1;
        for j = 1:4
            component = component*theta(k,j)^patterns(r,j)* ...
                (1-theta(k,j))^(1-patterns(r,j));
        end
        probabilities(r) = probabilities(r)+mixing(k)*component;
    end
end
probabilities = probabilities/sum(probabilities);
end
