

clear all, 
clc


%%%  Gloptipoly Code
mpol x 10

M = 500;

eta_true = [0.5, 0.5];

% pi_true(k, j, val+1) = P(y_j=val | class=k)
%          j=1: P(y=0), P(y=1) | j=2: P(y=0), P(y=1)
pi_true = cat(3, [0.4, 0.6; 0.8, 0.2], ...  
                 [0.1, 0.9; 0.6, 0.4]);

Ex_6_7_3_data                                              

N(1) = sum(Y(:, 1) == 0 & Y(:, 2) == 0); % [0 0] 
N(2) = sum(Y(:, 1) == 0 & Y(:, 2) == 1); % [0 1] 
N(3) = sum(Y(:, 1) == 1 & Y(:, 2) == 0); % [1 0] 
N(4) = sum(Y(:, 1) == 1 & Y(:, 2) == 1); % [1 1]
NN= N;
N = N/max(N);

p{1} = x(1)*x(3)*x(5) + x(2)*x(7)*x(9);
p{2} = x(1)*x(3)*x(6) + x(2)*x(7)*x(10);
p{3} = x(1)*x(4)*x(5) + x(2)*x(8)*x(9);
p{4} = x(1)*x(4)*x(6) + x(2)*x(8)*x(10);

Ktheta = [x(1)+x(2) == 1;
    x(3)+x(4) == 1;
    x(5)+x(6) == 1;
    x(7)+x(8) == 1;
    x(9)+x(10) == 1;
    x(1) >= 0; x(2) >= 0; x(3) >= 0; x(4) >= 0; x(5) >= 0; x(6) >= 0;
    x(7) >= 0; x(8) >= 0; x(9) >= 0; x(10) >= 0];

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

settings = sdpsettings('solver', 'mosek');
optimize(MomRelax,obj,settings); % Unreported warm-up solve
started = tic;
sol = optimize(MomRelax,obj,settings);
elapsed = toc(started)

mdim = MCone.f ;
for k=1:1
    locM = mat( c_Aty(mdim+1:mdim+(MCone.s(k))^2) );  
    for r = 1 : ord
        mlen_1 = nchoosek(length(x)+r,r);  
        eigs(value( locM(1:mlen_1,1:mlen_1) ), 5)',
    end
end 

moment = value(moment);

pp1 = value( b0{1}+bp{1}'*mom_y ),
pp2 = value( b0{2}+bp{2}'*mom_y ),
pp3 = value( b0{3}+bp{3}'*mom_y ),
pp4 = value( b0{4}+bp{4}'*mom_y ),


upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3) + NN(4)*log(pp4)

