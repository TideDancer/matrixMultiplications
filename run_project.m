% random matrix test 
A_norm = [];
B_norm = [];
C_norm = [];

range = 9:13;

for i = range
  i
  dim = 2^i;
  
  A = squareMatrixGen(dim, 'dense', 'normal');
  B = squareMatrixGen(dim, 'dense', 'normal');
  
  delta = 1/dim;
  epsilon = 1e-1;
  const = 1; 
  
  tic;
  C_approx1 = randomProjMult(A, B, 'kyrillidis2014approximate', [delta, epsilon, const]); % kyrillidis2014 paper
  toc;
  tic;
  C_approx2 = randomProjMult(A, B, 'FJLT', [delta, epsilon, const]);
  toc;
  tic;
  C_approx3 = randomProjMult(A, B, 'tug-of-war', [delta, epsilon, const]);
  toc;

  % ------------------- compare --------------------------
end


