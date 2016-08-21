% define dimension when run
dim = 2^9;

delta = 1/dim;  % failure probability
epsilon = 1e-1; % error norm <= epsilon ||A|| ||B||
beta = 1;
const = 1;

density = 'dense';
distribution = 'normal';
matrix = 'squareMatrixGen';

A = squareMatrixGen(dim, 'dense', 'normal');
B = squareMatrixGen(dim, 'dense', 'normal');
tic;
C = A*B;
toc;
nnzAB = nnz(C);
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');

