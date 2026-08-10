function out = Table_5_6_solve_log_relaxation(instance, use_lme, ord)
%TABLE_6_7_SOLVE_LOG_RELAXATION Solve a standard or LME log relaxation.
% This helper is called by Table_5.m and Table_6.m. Do not run it
% directly. It constructs the selected relaxation and delegates the
% moment-SDP calculation to Common_solve_log_moment_relaxation.m.
% Solve the standard or LME log-moment relaxation for a generated instance.
%
% J. Choi, August 4, 2026

mset clear
mset('verbose', false)
mpol('x', instance.n)

p = cell(instance.m, 1);
for j = 1:instance.m
    affine = instance.beta(j);
    for i = 1:instance.n
        affine = affine+instance.U(j, i)*x(i);
    end
    p{j} = instance.epsilon+affine^2;
end

Ktheta = [x'*x <= 1; 1-sum(x) >= 0];
for i = 1:instance.n
    Ktheta = [Ktheta; x(i) >= 0]; %#ok<AGROW>
end

flat_degree = 1;
if use_lme
    flat_degree = instance.m+1;
    g = cell(instance.n, 1);
    for i = 1:instance.n
        g{i} = 0;
        for j = 1:instance.m
            product_others = 1;
            for ell = 1:instance.m
                if ell ~= j
                    product_others = product_others*p{ell};
                end
            end
            g{i} = g{i}+instance.weights(j)*diff(p{j}, x(i))*product_others;
        end
    end
    % g{i} is p^1 times the ith partial derivative of the log objective.
    hat_tau0 = 0;
    for i = 1:instance.n
        hat_tau0 = hat_tau0+x(i)*g{i};
    end
    Ktheta = [Ktheta; hat_tau0 >= 0; ...
        hat_tau0*(1-sum(x)) == 0];
    for i = 1:instance.n
        hat_tau_i = hat_tau0-g{i};
        Ktheta = [Ktheta; hat_tau_i >= 0; ...
            hat_tau_i*x(i) == 0]; %#ok<AGROW>
    end
end

out = Common_solve_log_moment_relaxation(x, p, instance.weights, instance.weights, ...
    Ktheta, ord, ord, flat_degree, flat_degree);
end
