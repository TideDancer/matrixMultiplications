% random matrix test 
A_norm = [];
B_norm = [];
C_norm = [];
errRatio1 = []; errRatio2 = []; errRatio3 = [];

range = 12:12;

for i = range
  i
  dim = 2^i;
  
  A = squareMatrixGen(dim, 'dense', 'normal');
  B = squareMatrixGen(dim, 'dense', 'normal');
  
  delta = 1/dim;
  epsilon = 1e-1;
  const = 1; 
  
%  tic;
%  C_approx1 = randomProjMult(A, B, 'kyrillidis2014approximate', [delta, epsilon, const]); % kyrillidis2014 paper
%  toc;
  tic;
  C_approx1 = randomProjMult(A, B, 'clarkson2009numerical', [delta, epsilon, const]); % clarkson2009numerical paper
  toc;
 
  tic;
%  C_approx2 = randomProjMult(A, B, 'FJLT', [delta, epsilon, const]);
  toc;
  tic;
%  C_approx3 = randomProjMult(A, B, 'tug-of-war', [delta, epsilon, const]);
  toc;

  % ------------------- compare --------------------------
  tic;
  C = A*B;
  toc;
  A_norm = norm(A, 'fro');
  B_norm = norm(B, 'fro');
  AB_norm = A_norm * B_norm;
  C_norm = norm(C, 'fro');
  
  error1 = C - C_approx1;
  error1_norm = norm(error1, 'fro');
  errRatio1 = [errRatio1 error1_norm/AB_norm];
  
%  error2 = C - C_approx2;
%  error2_norm = norm(error2, 'fro');
%  errRatio2 = [errRatio2 error2_norm/AB_norm];
%
 % error3 = C - C_approx3;
 % error3_norm = norm(error3, 'fro');
 % errRatio3 = [errRatio3 error3_norm/AB_norm];
%%
end


