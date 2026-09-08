
clear all, 
clc

%%%  Gloptipoly Code
mpol x 3
p{1} = x(1)^2 + 2*x(1)*x(3);
p{2} = x(2)^2 + 2*x(2)*x(3);
p{3} = 2*x(1)*x(2);
p{4} = x(3)^2;

N(1) = 182; N(2) = 60; N(3) = 17; N(4) = 176;
NN= N;
N = N/max(N);

% Lagrange Multiplier (Simplex)
lme0 = 2*sum(N);
lme(1) = lme0*x(1)*(x(1)+2*x(3))-(N(1)+N(3))*(x(1)+2*x(3))-N(1)*x(1);
lme(2) = lme0*x(2)*(x(2)+2*x(3))-(N(2)+N(3))*(x(2)+2*x(3))-N(2)*x(2);
lme(3) = lme0*x(3)*(x(1)+2*x(3))*(x(2)+2*x(3))-2*N(4)*(x(1)+2*x(3))*(x(2)+2*x(3))-2*N(1)*x(3)*(x(2)+2*x(3))-2*N(2)*x(3)*(x(1)+2*x(3));

Ktheta = [1 >= sum(x)];

% Feasibility
Ktheta = [Ktheta; 1-x'*x >= 0];
for i = 1:length(x)
    Ktheta = [Ktheta; x(i) >= 0];
end

% LME Positivity
for i = 1:length(x)
    Ktheta = [Ktheta; lme(i) >= 0];
end

% Complementarity
Ktheta = [Ktheta; 1 - sum(x) == 0];
for i = 1:length(x)
    Ktheta = [Ktheta; lme(i)*x(i) == 0];
end

ord = 3;

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

settings = sdpsettings('solver', 'mosek');
warmup = optimize(MomRelax,obj,settings);
tic
sol = optimize(MomRelax,obj,settings);
solve_time = toc;

mdim = MCone.f + MCone.l;
for k=1:1
    locM = mat( c_Aty(mdim+1:mdim+(MCone.s(k))^2) );  
    for r = 1 : ord
        mlen_1 = nchoosek(length(x)+r,r);  
        eigs(value( locM(1:mlen_1,1:mlen_1) ), 10)',
    end
end 

moment = value(moment);

RR = rank(moment, 10e-5)

%uu = value( mom_y(1:4)' ),

pp1 = value( b0{1}+bp{1}'*mom_y ),
pp2 = value( b0{2}+bp{2}'*mom_y ),
pp3 = value( b0{3}+bp{3}'*mom_y ),
pp4 = value( b0{4}+bp{4}'*mom_y ),

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3) + NN(4)*log(pp4) ;
solve_time

%%% Extract supports
xx = extractmin(moment, length(x), ord);
xx{1}
