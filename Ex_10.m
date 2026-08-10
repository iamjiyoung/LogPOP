% Example 10.
%
% J. Choi, August 10, 2026

clear all, 
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))


%%%  Gloptipoly Code
mpol x 12

n = 12;

v = zeros(n, 10);
v(1:2, 1) = [-0.2; -1.6];
v(3:4, 2) = [-1; -0.4];
v(5:6, 3) = [0.4; -1];
v(7:8, 4) = [0.4; 0.8];
v(9:10, 5) = [-1.1; -1.5];
v(11:12, 6) = [-0.8; -0.7];
v([1, 7], 7) = [-0.3; 0];
v([2, 8], 8) = [-1.7; -1];
v([3, 9], 9) = [1.2; -1.9];
v([4, 10], 10)= [1.7; 0.9];

p{1} = 20 - 1.5*( (x(2)+v(1,1)) + (x(3)+v(2,1)) )^4;
p{2} = 25 - 1.2*( (x(4)+v(3,2)) + (x(5)+v(4,2)) )^4;
p{3} = 18 - 2.5*( (x(6)+v(5,3)) + (x(7)+v(6,3)) )^4;
p{4} = 22 - 1.1*( (x(8)+v(7,4)) + (x(9)+v(8,4)) )^4;
p{5} = 30 - 1.6*( (x(10)+v(9,5)) + (x(11)+v(10,5)) )^4;
p{6} = 15 - 1.9*( (x(1)+v(11,6)) - (x(6)+v(12,6)) + (x(11)-1.3) )^4;
p{7} = 28 - 2.3*( -(x(2)+v(1,7)) + (x(7)+v(7,7)) - (x(12)+0.7) )^4;
p{8} = 19 - 1.7*( (x(3)+v(2,8)) - (x(8)+v(8,8)) + (x(1)-0.7) )^4;
p{9} = 21 - 2.4*( -(x(4)+v(3,9)) + (x(9)+v(9,9)) - (x(2)+0.8) )^4;
p{10} = 26 - 1.4*( (x(5)+v(4,10)) - (x(10)+v(10,10)) + (x(3)-0.5) )^4;

% Objective weights
N = [1 1 1 1 1 1 1 1 1 1];
NN = N;
N = N/max(N);

% feasibility
Ktheta = [(sum(x(1:6)))^2 - 15 <= 0;
    (sum(x(7:12)))^2 - 15 <= 0;
    (x(1)-x(3)+x(5)-x(7)+x(9)-x(11))^2 - 8 <= 0;
    (x(2)-x(4)+x(6)-x(8)+x(10)-x(12))^2 - 8 <= 0;
    (x(1)+x(12))^2 + (x(2)+x(11))^2 - 9 <= 0;
    (x(3)-x(10))^2 + (x(4)-x(9))^2 + (x(5)-x(8))^2 - 9 <= 0;
    -sum(x) <= 0; % sum(x) >= 0
    sum(x(1:2:11)) - sum(x(2:2:12)) - 5 <= 0; % (x1-x2)+(x3-x4)+... <= 5
    (2*x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 + 3*x(5)^2 + x(6)^2 + ...
     2*x(7)^2 + x(8)^2 + 2*x(9)^2 + x(10)^2 + 3*x(11)^2 + x(12)^2) - 20 <= 0;
    (x(1)+x(6)+x(12))^2 - 4 <= 0;
    x(1) - x(12) - 3 <= 0;
    sum(x) >= 0;
    p{1} >= 0;
    p{2} >= 0;
    p{3} >= 0;
    p{4} >= 0;
    p{5} >= 0;
    p{6} >= 0;
    p{7} >= 0;
    p{8} >= 0;
    p{9} >= 0;
    p{10} >= 0; ];

ord = 2;

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
sol = optimize(MomRelax,obj, ...
    sdpsettings('solver', 'mosek', 'verbose', 0));
solve_time = toc(started);

mdim = MCone.f + MCone.l;
for k=1:1
    locM = mat( c_Aty(mdim+1:mdim+(MCone.s(k))^2) );  
    for r = 1 : ord
        mlen_1 = nchoosek(length(x)+r,r);  
        eigs(value( locM(1:mlen_1,1:mlen_1) ), 5);
    end
end 

moment = value(moment);

pp1 = value( b0{1}+bp{1}'*mom_y );
pp2 = value( b0{2}+bp{2}'*mom_y );
pp3 = value( b0{3}+bp{3}'*mom_y );
pp4 = value( b0{4}+bp{4}'*mom_y );
pp5 = value( b0{5}+bp{5}'*mom_y );
pp6 = value( b0{6}+bp{6}'*mom_y );
pp7 = value( b0{7}+bp{7}'*mom_y );
pp8 = value( b0{8}+bp{8}'*mom_y );
pp9 = value( b0{9}+bp{9}'*mom_y );
pp10 = value( b0{10}+bp{10}'*mom_y );

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3) + NN(4)*log(pp4) + NN(5)*log(pp5);
upperbound = upperbound + NN(6)*log(pp6) + NN(7)*log(pp7) ...
    + NN(8)*log(pp8) + NN(9)*log(pp9) + NN(10)*log(pp10);

ext = value(mom_y(1:12));
first_moment_value = NN(1)*log(double(subs(p{1},x,ext))) ...
    + NN(2)*log(double(subs(p{2},x,ext))) ...
    + NN(3)*log(double(subs(p{3},x,ext))) ...
    + NN(4)*log(double(subs(p{4},x,ext))) ...
    + NN(5)*log(double(subs(p{5},x,ext))) ...
    + NN(6)*log(double(subs(p{6},x,ext))) ...
    + NN(7)*log(double(subs(p{7},x,ext))) ...
    + NN(8)*log(double(subs(p{8},x,ext))) ...
    + NN(9)*log(double(subs(p{9},x,ext))) ...
    + NN(10)*log(double(subs(p{10},x,ext)));

[xx, extract_t, extract_rank] = extractmin_flat( ...
    moment, length(x), ord, 1e-4, 2, true, 2);

fprintf('\n=== Example 10: reported result ===\n');
fprintf('Relaxation order k        : %d\n', ord);
fprintf('Optimal relaxation value  : %.10f\n', upperbound);
fprintf('Flat-truncation order t    : %d\n', extract_t);
fprintf('Rank                       : %d\n', extract_rank);
fprintf('Runtime (seconds)          : %.4f\n', solve_time);
fprintf('Extracted optimizer        : %s\n', mat2str(xx{1}(:)', 10));

%%% Verify the strict-positivity witness used in the paper.
witness = [-0.3; -0.9; 1.5; -0.9; 1.1; -0.4; ...
    0.5; -1.1; 1.0; 2.1; 0.6; 0.7];
witness_p = zeros(1, length(p));
for i = 1:length(p)
    witness_p(i) = double(subs(p{i}, x, witness));
end

witness_slacks = [
    15-sum(witness(1:6))^2;
    15-sum(witness(7:12))^2;
    8-(witness(1)-witness(3)+witness(5)-witness(7)+witness(9)-witness(11))^2;
    8-(witness(2)-witness(4)+witness(6)-witness(8)+witness(10)-witness(12))^2;
    9-(witness(1)+witness(12))^2-(witness(2)+witness(11))^2;
    9-(witness(3)-witness(10))^2-(witness(4)-witness(9))^2-(witness(5)-witness(8))^2;
    sum(witness);
    5-(sum(witness(1:2:11))-sum(witness(2:2:12)));
    20-(2*witness(1)^2+witness(2)^2+2*witness(3)^2+witness(4)^2+ ...
        3*witness(5)^2+witness(6)^2+2*witness(7)^2+witness(8)^2+ ...
        2*witness(9)^2+witness(10)^2+3*witness(11)^2+witness(12)^2);
    4-(witness(1)+witness(6)+witness(12))^2;
    3-witness(1)+witness(12)];

assert(min(witness_p) > 0, ...
    'The paper witness does not make every logarithm argument positive.');
assert(min(witness_slacks) >= -1e-10, ...
    'The paper witness violates a defining constraint.');

[minimum_p, minimum_p_index] = min(witness_p);
fprintf('\n=== Example 10: positivity-witness check ===\n');
fprintf('Minimum p_i at witness     : %.10f (i = %d)\n', ...
    minimum_p, minimum_p_index);
fprintf('Minimum constraint slack   : %.10f\n', min(witness_slacks));
