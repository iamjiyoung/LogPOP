% LME relaxation in Example 6(i).
%
% J. Choi, July 31, 2026

clear all,
clc
run(fullfile(fileparts(mfilename('fullpath')), ...
    'functions', 'setup_paths.m'))

%%%  Gloptipoly Code
mpol x 3
p{1} = x(1)^3 + 3*x(1)^2*x(2) + 3*x(1)^2*x(3) ;
p{2} = 3*x(1)*x(2)^2 + 6*x(1)*x(2)*x(3);
p{3} = 3*x(1)*x(3)^2;
p{4} = x(2)^3 + 3*x(2)^2*x(3);
p{5} = 3*x(2)*x(3)^2;
p{6} = x(3)^3;

%Choose random N
%ord 3
N=[0.0968    0.1419    0.2194    0.0839    0.2839    0.1742]; 

%ord 3
% N=[0.7691    0.6389    0.8931    0.0607    0.1758    0.4163];

% %ord 3
% N=[0.1774    0.3959    0.4922    0.4379    0.6354    0.1527]; 
 
% %ord 3
% N=[   0.2920    0.4317    0.0155    0.9841    0.1672    0.1062];
% 
% %ord 3
% N=[  0.0835    0.6260    0.6609    0.7298    0.8908    0.9823];
% 
% %ord 3
% N=[0.3424    0.7360    0.7947    0.5449    0.6862    0.8936];


% Calculate gradient (here g{i} = g{i}*x(i))
g{1} = N(1)*3*(x(1)+2*x(2)+2*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) + (N(2) + N(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));
g{2} = N(1)*3*x(2)*(x(2)+2*x(3))*(x(2)+3*x(3)) + N(2)*(2*x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3)) + N(4)*3*(x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3)) + N(5)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));
g{3} = N(1)*3*x(3)*(x(2)+2*x(3))*(x(2)+3*x(3)) + N(2)*2*x(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3)) + 2*N(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) + N(4)*3*x(3)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3)) + (2*N(5) +3*N(6))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3));

% Lagrange Multiplier (Simplex)
lme0 = g{1} + g{2} + g{3};
lme(1) = lme0*x(1) - (N(1)*3*(x(1)+2*x(2)+2*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) + (N(2) + N(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)));
lme(2) = lme0*x(2) - (3*N(1)*(x(2)+2*x(3))*(x(2)+3*x(3))*x(2) + N(2)*(2*x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3)) + N(4)*3*(x(2)+2*x(3))*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3)) + N(5)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)));
lme(3) = lme0*x(3) - (3*N(1)*(x(2)+2*x(3))*(x(2)+3*x(3))*x(3) + N(2)*2*(x(1)+3*x(2)+3*x(3))*(x(2)+3*x(3))*x(3) + N(3)*2*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) + N(4)*3*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*x(3) + N(5)*2*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)) + 3*N(6)*(x(1)+3*x(2)+3*x(3))*(x(2)+2*x(3))*(x(2)+3*x(3)));

% Feasibility
Ktheta = [x'*x <= 1;
    1-sum(x) >= 0];
for i = 1:length(x)
    Ktheta = [Ktheta; x(i) >= 0];
end

% LME Positivity
Ktheta = [Ktheta; lme0 >= 0];
for i = 1:length(x)
    Ktheta = [Ktheta; lme(i) >= 0];
end

% Complementarity
for i = 1:length(x)
    Ktheta = [Ktheta; lme(i)*x(i) == 0];
end
Ktheta = [Ktheta; lme0*(1-sum(x)) == 0];

% d = 3;
ord = 4;

for i = 1:length(N)
    PO{i} = msdp (max(p{i}), Ktheta, ord);
    [A{i}, bp{i}, c{i}, K{i}, b0{i}, s{i}] = msedumi(PO{i});
end

%%%% YALMIP Code

MCone = K{1};
mom_cc = c{1};
len_y = length( bp{1} ); 
mom_y = sdpvar( len_y, 1);
Aty = A{1}'*mom_y; 
c_Aty = c{1}-Aty; 
MomRelax = []; 
for k=1:MCone.f
    MomRelax = [MomRelax, c_Aty(k)== 0 ];
end
mdim = MCone.f ;
for k=1:MCone.l
    MomRelax = [MomRelax, c_Aty(mdim+k)>=0];
end
mdim = mdim + MCone.l;
for k=1:length( MCone.s )
    locM = mat( c_Aty(mdim+1:mdim+(MCone.s(k))^2) );  
    if k == 1
        moment = locM;
    end
    MomRelax = [MomRelax, locM >= 0 ];
    mdim = mdim + ( MCone.s(k) )^2;
end 

obj = - N(1)*log(b0{1}+bp{1}'*mom_y);
for i = 2:length(p)
    obj = obj - N(i)*log(b0{i}+bp{i}'*mom_y);
end

started = tic;
sol = optimize(MomRelax,obj, sdpsettings('solver', 'mosek'));
solve_time = toc(started);

mdim = MCone.f + MCone.l;
for k=1:1
    locM = mat( c_Aty(mdim+1:mdim+(MCone.s(k))^2) );  
    for r = 1 : ord
        mlen_1 = nchoosek(length(x)+r,r);  
        eigs(value( locM(1:mlen_1,1:mlen_1) ), 5);
    end
end 
%%
mlen_1 = nchoosek(length(x)+ord,ord); 
moment = value( locM(1:mlen_1,1:mlen_1) );

pp1 = value( b0{1}+bp{1}'*mom_y );
pp2 = value( b0{2}+bp{2}'*mom_y );
pp3 = value( b0{3}+bp{3}'*mom_y );
pp4 = value( b0{4}+bp{4}'*mom_y );
pp5 = value( b0{5}+bp{5}'*mom_y );
pp6 = value( b0{6}+bp{6}'*mom_y );

upperbound = N(1)*log(pp1) + N(2)*log(pp2) + N(3)*log(pp3) ...
    + N(4)*log(pp4) + N(5)*log(pp5) + N(6)*log(pp6);
%%
%%% Extract supports
[xx, extract_t, extract_rank] = extractmin_flat( ...
    moment, length(x), ord, 1e-4, 3, true, 3);

fprintf('\n=== Example 6(i): LME moment relaxation ===\n');
fprintf('Relaxation order k        : %d\n', ord);
fprintf('Optimal relaxation value  : %.10f\n', upperbound);
fprintf('Flat-truncation order t    : %d\n', extract_t);
fprintf('Rank                       : %d\n', extract_rank);
fprintf('Runtime (seconds)          : %.4f\n', solve_time);
fprintf('Extracted atoms:\n');
for atom_index = 1:numel(xx)
    fprintf('  atom %d: %s\n', atom_index, mat2str(xx{atom_index}(:)', 10));
    p_at_atom = zeros(1, numel(p));
    for p_index = 1:numel(p)
        p_at_atom(p_index) = double(subs(p{p_index}, x, xx{atom_index}));
    end
    objective_at_atom = sum(N.*log(p_at_atom));
    fprintf('  objective at atom %d     : %.10f\n', ...
        atom_index, objective_at_atom);
    fprintf('  min_i p_i(atom %d)       : %.12e\n', ...
        atom_index, min(p_at_atom));
end
