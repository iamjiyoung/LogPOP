function [x, t, r] = extractmin_flat( ...
        M, nbvar, maxdeg, tol, mindeg, warn_missing, degree_lag)
% Extract atoms when rank(M_t) = rank(M_{t-degree_lag}).

if nargin < 4
    tol = 1e-4;
end
if nargin < 5
    mindeg = 1;
end
if nargin < 6
    warn_missing = true;
end
if nargin < 7
    degree_lag = 1;
end
if degree_lag < 1 || degree_lag ~= floor(degree_lag)
    error('extractmin_flat:InvalidDegreeLag', ...
        'degree_lag must be a positive integer.');
end

M = (M + M')/2;
expected = nchoosek(nbvar + maxdeg, maxdeg);
if size(M,1) < expected || size(M,2) < expected
    error('extractmin_flat:InvalidSize', ...
        'Moment matrix is too small for %d variables and order %d.', ...
        nbvar, maxdeg);
end

x = {};
t = NaN;
r = NaN;
for candidate = maxdeg:-1:max(mindeg, degree_lag)
    ncurr = nchoosek(nbvar + candidate, candidate);
    previous_order = candidate-degree_lag;
    nprev = nchoosek(nbvar + previous_order, previous_order);
    Mcurr = M(1:ncurr, 1:ncurr);
    rcurr = rank(Mcurr, tol);
    rprev = rank(M(1:nprev, 1:nprev), tol);
    if rcurr == rprev
        x = extractmin(Mcurr, nbvar, candidate);
        t = candidate;
        r = rcurr;
        return
    end
end

if warn_missing
    warning('extractmin_flat:NoFlatTruncation', ...
        ['No flat truncation satisfying rank(M_t) = ' ...
        'rank(M_{t-%d}) was found up to order %d.'], ...
        degree_lag, maxdeg);
end
end
