

clear all, 
clc

a = [0;0;0;0;0];
b = 1;

%%%  Gloptipoly Code
mpol x 5

p_11 = x(3)+x(4);
p_12 = 1+2*x(1)+3*x(2)+3*x(4)+2*x(5);
p_21 = x(2)+x(4)-x(5);
p_22 = 2*x(2)+3*x(5);
p_31 = x(1)+x(3)+x(4);
p_32 = x(1)-x(3)+x(4);

p{1} = (p_11)^2 + (p_12)^2 + 0.01;
p{2} = (p_21)^2 + (p_22)^2 + 0.02;
p{3} = (p_31)^2 + (p_32)^2 + 0.03;


%Choose random N
N = [30/218    97/218    91/218];
NN = N;
N = N/max(N);

% feasibility
Ktheta = [b^2 - (x-a)'*(x-a) >= 0];

% d = 1;
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
warmup = optimize(MomRelax,obj,settings);
tic
sol = optimize(MomRelax,obj,settings);
toc

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

upperbound = NN(1)*log(pp1) + NN(2)*log(pp2) + NN(3)*log(pp3) ;

%%% Extract supports
xx = extractmin(moment, length(x), ord);
xx{1}

