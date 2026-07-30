function  x = extractmin(M,nbvar,dgr)

    
    % *****************************************
    % Compute SVD of moment matrices
    % to detect global optimality

    [U, Lmd]=eig(M);
    
    Lmd = diag(Lmd);
    n = length(Lmd);
    drop = find(Lmd>=1e-4);
    %drop = find(Lmd>=1e-6);
    
    rkM = length(drop);
    
    % solutions
    x = {};
    
    % ********************************************************************
    % Extract solutions
    %
    % Algorithm:
    % 1. extract a Cholesky factor U of moment matrix SED.M{K} such that
    %    U*U' = SED.M{K}
    % 2. reduce U to column echelon form via Gaussian elimination with
    %    column pivoting, identify monomial basis B(1)..B(P) in U (pivots)
    % 3. for each monomial X(I), I=1..P extract the coefficient
    %    matrix N(I) of monomials X(I)*B(1)..X(I)*B(P) in basis B(1)..B(P)
    % 4. compute common roots of multiplication matrices N(1)..N(P) and store
    %    them in output cell array X
    
    

    % Step 1: extract Cholesky factor of moment matrix
      
    U = U(:,drop)*diag(sqrt(Lmd(drop)));
    %norm(M-U*U');
     
    
    % Step 2: reduce Cholesky factor to column echelon form
    % and identify monomial basis
    [U,basis] = cef(U,1e-6);
   

    % Step 3: extract multiplication matrix for each variabl
    N = cell(1,nbvar);
    mon = mon_pow_vec(nbvar,dgr); % powers of monomials
    nmon = size(mon,1); % total number of monomials
    
    Bpow =[];
    for j = 1: nmon
      od = mon_order(mon(j,:));
      inbasis = find( basis == od);
      if ~isempty(inbasis), Bpow=[Bpow; mon(j,:)]; end
    end
    
    
    for i =1 : nbvar
    ei = eye(nbvar);
    ei = ei(i,:);
       row=[];
       for j = 1: size(Bpow,1);
        row = [row mon_order(Bpow(j,:)+ei)];
       end
       row = sort(row);
    N{i}=U(row,:);
    end
 
    
    % random combinations of multiplication matrices
      coef = rand(nbvar,1); coef = coef / sum(coef);
      MM = zeros(rkM);
      for i = 1:nbvar,
	    MM = MM + coef(i)*N{i};
      end;
	
      % ordered Schur decomposition of M      
      [Q,T] = orderschur(MM);
	
      % retrieve optimal vectors
      % it is assumed than there is no multiple root
      x = cell(1,rkM);
    
      for i = 1:rkM,
	   x{i} = zeros(nbvar,1);
	    for j = 1:nbvar,
	    x{i}(j) = Q(:,i)'*N{j}*Q(:,i);
	    end;
      end
      
      %x{:} 
      

return;

% ********************************************************************
% Utilities

function mat = filtermat(mat,typevar,base)
% Filter 0-1 and +/-1 constraints in matrix MAT
% according to type variable vector TYPEVAR and base BASE
if any(typevar),
 nbvar = length(typevar);
 % Decomposition of MAT in base BASE
 mask(:,:,nbvar) = mat;
 for k = nbvar:-1:2,
   remainder = rem(mask(:,:,k), base);
   mask(:,:,k-1) = (mask(:,:,k)-remainder)/base;
   mask(:,:,k) = remainder;
 end;
 % Recomposition with filtering
 mat = zeros(size(mat));
 for k = 1:nbvar,
  mat = mat * base;
  switch typevar(k),
   case 1, % 0-1 constraint : x^2 -> x
    mat = mat + (mask(:,:,k) ~= 0); 
   case -1, % +/-1 constraint: x^2 -> 0
    mat = mat + rem(mask(:,:,k), 2);
   case 0, % no constraint
    mat = mat + mask(:,:,k);
   otherwise
    error('Invalid entry in variable type vector');
  end;
 end;
end;

function multi = locate(indices,size)
% Transform linear indices INDICES into a cell array of
% multidimensional indices in a matrix of size SIZE
% For example, LOCATE(5,[2 3]) returns [1 3], i.e. 5th element in
% 2x3 matrix is located in first row and third column 
if any(size),
 multi = cell(length(indices),1);
 for k1 = 1:length(indices),
  index = indices(k1);
  for k2 = 1:length(size),
   multi{k1}(k2) =  1+rem(index-1,size(k2));
   index = 1+floor((index-1)/size(k2));
  end;
 end;
else
 multi = {1};
end;

function t = generate(digits,sum,base)
% Generate all the indices with DIGITS digits summing to SUM in base BASE.
% For example, GENERATE(3,2,3) returns the vector [2 4 10 6 12 18] which is
% [002 011 101 020 110 200] in base 3. Note that the most significant
% digit is located at the left.  
if digits < 1,
  t = [];
elseif digits < 2,
  t = sum;
else,
  t = zeros(1,nchoosek(sum+digits-1,digits-1));
  j = 1;
  for i = sum:-1:0,
    s = generate(digits-1,sum-i,base); r = length(s); % recursive call
    t(j:j+r-1) = i*ones(1,r)+s*base; j = j+r;
  end;
end;

function [A,basis] = cef(A,tol)
% The instruction
%
%  [E,BASIS] = CEF(A)
%
% computes a column echelon form E of matrix A
% and returns basis row indices in vector BASIS
%
% The relative threshold for column pivoting can be specified
% as an additional input argument
%
% The reduction is performed by Gaussian elimination with
% column pivoting, based on Matlab's RREF routine
  
[n,m] = size(A);

% Loop over the entire matrix.
i = 1; j = 1; basis = [];
while (i <= m) & (j <= n)
   % Find value and index of largest element in the remainder of row j
   [p,k] = max(abs(A(j,i:m))); k = k+i-1;
   if (p <= tol)
      % The row is negligible, zero it out.
      A(j,i:m) = zeros(1,m-i+1,1);
      j = j + 1;
   else
      % Remember row index
      basis = [basis j];
      % Swap i-th and k-th columns
      A(j:n,[i k]) = A(j:n,[k i]);
      % Find a non-negligible pivot element in the column
      found = 0;
      while ~found,
	if abs(A(j,i)) < tol*max(abs(A(:,i)))
	  j = j + 1;
	  found = (j == n);
	else
	  found = 1;
	end;
      end;
      if j <= n,
        % Divide the pivot column by the pivot element
	A(j:n,i) = A(j:n,i)/A(j,i);
	% Subtract multiples of the pivot column from all the other columns
	for k = [1:i-1 i+1:m]
	  A(j:n,k) = A(j:n,k) - A(j,k)*A(j:n,i);
	end
	i = i + 1;
	j = j + 1;
      end
   end
end

function P = evaluate(P,x)
% The instruction VAL = EVALUATE(P,X) evaluates an objective function,
% inequality or equality constraint P (in GloptiPoly's format) at a
% vector X

if issparse(P.c) % sparse coefficient matrix

 nd = length(P.s);
 d = P.s;
 val = 0;
 nzndx = find(P.c);
 for j = 1:length(nzndx),
  k = [1 cumprod(d(1:nd-1))];
  ndx = nzndx(j) - 1;
  pow = zeros(1,nd);
  for i = nd:-1:1,
   pow(i) = floor(ndx/k(i));
   ndx = rem(ndx,k(i));
  end;
  coef = 1;
  for i = 1:nd,
   coef = coef*x(i)^pow(i);
  end;
  val = val + P.c(nzndx(j))*coef;
 end;
 P = val;
 
else % non-sparse matrix, applies Horner's scheme

 P = P.c;
 nd = ndims(P);
 if (nd == 2) & (size(P,2) == 1), nd = 1; end;
 d = size(P);
 for n = 1:nd, ind{n} = 1:d(n); end;
 for n = nd:-1:1,
  indP = {ind{1:n-1} d(n)};
  newP = P(indP{:});
  for k = d(n)-1:-1:1,
   indP = {ind{1:n-1} k};
   newP = P(indP{:}) + x(n)*newP;
  end;
  P = newP;
 end;

end;

function [U,T] = orderschur(X)
% [U,T] = ORDERSCHUR(X) computes the real Schur decomposition X = U*T*U'
% of a matrix X with real eigenvalues sorted in increasing order along
% the diagonal of T
%
% Algorithm: perform unordered Schur decomposition with Matlab's SCHUR
% function, then order the eigenvalues with Givens rotations, see
% [Golub, Van Loan. Matrix computations. 1996]

[U,T] = schur(X);
U = real(U); T = real(T);
n = size(X,1);

order = 0;
while ~order, % while the order is not correct
 order = 1;
 for k = 1:n-1,
  if T(k,k) - T(k+1,k+1) > 0,
   order = 0; % diagonal elements to swap
   % Givens rotation
   [c,s] = givens(T(k,k+1),T(k+1,k+1)-T(k,k));
   T(k:k+1,k:n) = [c s;-s c]'*T(k:k+1,k:n);
   T(1:k+1,k:k+1) = T(1:k+1,k:k+1)*[c s;-s c];
   U(1:n,k:k+1) = U(1:n,k:k+1)*[c s;-s c];
  end;
 end;
end; % while

function [c,s] = givens(a,b)
% Givens rotation for ordered Schur decomposition
if b == 0
  c = 1; s = 0;
else
  if abs(b) > abs(a)
    t = -a/b; s = 1/sqrt(1+t^2); c = s*t;
  else
    t = -b/a; c = 1/sqrt(1+t^2); s = c*t;
  end
end






