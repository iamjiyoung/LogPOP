% Example 9.
%
% J. Choi, July 29, 2026

clear all, 
clc
run(fullfile(fileparts(mfilename('fullpath')), 'functions', 'setup_paths.m'))


%%%  Gloptipoly Code
mpol x 10

p{1} = 10 - ( 2*(x(1)+0.5)^2 + 3*(x(2)+0.5)^2 -3* (x(1)-0.5)*(x(2)+0.5) );
p{2} = 12 - ( 3*(x(3)-0.6)^2 + 2*(x(4)+0.1)^2 - (x(3)-0.6)*(x(4)-0.1) );
p{3} = 11 - ( 4*(x(5)+0.7)^2 + (x(6)-0.2)^2 );
p{4} = 15 - ( 2*(x(7)-0.8)^2 + 5*(x(8)+0.3)^2 - 2*(x(7)-0.8)*(x(8)+0.3) );
p{5} = 13 - ( -3*(x(9)-0.9)^2 + 2*(x(10)+0.4)^2 );

%Choose random N
N = [3 2 1 1 5];
NN = N;
N = N/max(N);

% feasibility
Ktheta = [1 - p{1} <= 0;  
    2 - p{2} <= 0;  
    3 - p{3} <= 0;  
    4 - p{4} <= 0;  
    5 - p{5} <= 0; 
    (x(1)+x(2)+x(3))^2 + (x(4)+x(5))^2 - 8 >= 0;
    (x(6)-x(7)-x(8))^2 + x(9)^2 - 6 <= 0;
    (x(1)-x(10))^2 + (x(2)-x(9))^2 - 7 <= 0;
    (x(2)+x(4)-x(6)+x(8)-x(10))^2 - 9 >= 0];

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

moment = value(moment);

pp1 = value( b0{1}+bp{1}'*mom_y );
pp2 = value( b0{2}+bp{2}'*mom_y );
pp3 = value( b0{3}+bp{3}'*mom_y );
pp4 = value( b0{4}+bp{4}'*mom_y );
pp5 = value( b0{5}+bp{5}'*mom_y );

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3) ...
    + NN(4)*log(pp4) + NN(5)*log(pp5);

[xx, extract_t, extract_rank] = extractmin_flat(moment, length(x), ord);

fprintf('\n=== Example 9: reported result ===\n');
fprintf('Relaxation order k        : %d\n', ord);
fprintf('Optimal relaxation value  : %.10f\n', upperbound);
fprintf('Flat-truncation order t    : %d\n', extract_t);
fprintf('Rank                       : %d\n', extract_rank);
fprintf('Runtime (seconds)          : %.4f\n', solve_time);
fprintf('Extracted optimizer        : %s\n', mat2str(xx{1}(:)', 10));
