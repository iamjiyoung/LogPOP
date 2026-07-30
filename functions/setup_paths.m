% J. Choi, July 29, 2026
%
% Add the numerical packages and local function directories to the MATLAB path.

functions_dir = fileparts(mfilename('fullpath'));
logpop_root = fileparts(functions_dir);
candidate_roots = {fileparts(logpop_root), ...
    fileparts(fileparts(logpop_root))};
packages = '';

for path_index = 1:numel(candidate_roots)
    candidate = fullfile(candidate_roots{path_index}, '00_Packages');
    if exist(candidate, 'dir')
        packages = candidate;
        break
    end
end

if isempty(packages)
    error(['The 00_Packages directory was not found beside real_codes ' ...
        'or one directory above it.']);
end

required = {
    fullfile(packages, 'gloptipoly3'), ...
    fullfile(packages, 'YALMIP-master'), ...
    fullfile(packages, 'SeDuMi_1_3'), ...
    fullfile(packages, 'Mosek', '10.0', 'toolbox', 'r2017a')};

for path_index = 1:numel(required)
    if ~exist(required{path_index}, 'dir')
        error('Required package directory not found: %s', ...
            required{path_index});
    end
end

addpath(genpath(required{1}));
addpath(genpath(required{2}));
addpath(required{3});
addpath(required{4});
addpath(functions_dir);
addpath(fullfile(functions_dir, 'extraction'));

clear functions_dir logpop_root candidate_roots candidate packages required path_index
