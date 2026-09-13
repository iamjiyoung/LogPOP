

clear all, 
clc


%%%  Gloptipoly Code
mpol x 6

M = 500;

eta_true = [0.5, 0.5];
pi_true = [0.4, 0.6;  % Class 1: P(y=0), P(y=1)
           0.8, 0.2]; % Class 2: P(y=0), P(y=1)

Ex_6_7_1_data

N(1) = 500-nnz(Y);
N(2) = nnz(Y);
NN= N;
N = N/max(N);

p{1} = x(1)*x(3) + x(2)*x(5);
p{2} = x(1)*x(4) + x(2)*x(6);

Ktheta = [x(1)+x(2) == 1;
    x(3)+x(4) == 1;
    x(5)+x(6) == 1;
    x(1) >= 0; x(2) >= 0; x(3) >= 0; x(4) >= 0; x(5) >= 0; x(6) >= 0];

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

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) 

%%% Extract supports
% xx = extractmin(moment, length(x), ord);
% xx{1}
% xx{2}

