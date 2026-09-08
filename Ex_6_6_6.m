
clear all, 
clc

%%%  Gloptipoly Code
mpol x 1
p{1} = 0.25*x(1) + 0.5*(1-x(1));
p{2} = 0.25*x(1) + 0.5*(1-x(1));
p{3} = 0.25*x(1);
p{4} = 0.25*x(1);


NN = [   59     8    16    17];
N = NN/max(NN); 

for i = 1:length(p)
    prod = 1;
    for j = 1:length(p)
        if j ~= i
            prod = prod*p{j};
        end
    end
    lcq{i} = prod;
end
% Calculate gradient 
for i = 1:length(x)
    summation = 0;
    for j = 1:length(N)
        summation = summation + N(j)*lcq{j}*diff(p{j}, x(i));
    end
    g{i} = summation;
end

% Lagrange Multiplier (Simplex)
lme0 = x(1)*g{1}  ;
lme = lme0 - [g{1}]';

Ktheta = [1-x'*x >= 0; 
    sum(x)-1 <= 0];

% Feasibility
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


RR = rank(moment, 10e-5)


pp1 = value( b0{1}+bp{1}'*mom_y ),
pp2 = value( b0{2}+bp{2}'*mom_y ),
pp3 = value( b0{3}+bp{3}'*mom_y ),
pp4 = value( b0{4}+bp{4}'*mom_y ),

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3)  + NN(4)*log(pp4) ;

%%% Extract supports
xx = extractmin(moment, length(x), ord);
xx{1}
1-xx{1}

