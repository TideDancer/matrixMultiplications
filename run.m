% ############# we can use ||C - C_approx || as a measurement of elementwise error #############################

% random matrix test 
dim = 2^10;

A = squareMatrixGen(dim, 'dense', 'normal');
B = squareMatrixGen(dim, 'dense', 'normal');

C_approx = basicMatrixMult(A, B, 'column2norm', 1); % coloum 2-norm based sampling

C_approx = elementMatrixMult(A, B, 'l2', 1); % l-2 based element-wise sampling

C_approx = randomProjMult(A, B, 'kyrillidis2014approximate', [1e-3, 1e-3, 1e-6]); % kyrillidis2014 paper
C_approx = randomProjMult(A, B, 'FJLT', [1e-2, 1e-2, 1]);

C_approx = nystrom(A, B);

C_approx = compressedSensing(A, B, [1, 1]);

C_approx = equationSolve(A, B);

C_approx = frequencyCounting(A, B); % A, B need to be non-negative matrix !!!!!!!!!!

C_approx = compressedFFT(A, B, 20); % parameter default value is 20




% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');

