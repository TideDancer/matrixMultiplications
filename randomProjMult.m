% matrix multiplication using random projection
% input A m*n, B n*p
m = 1024;
n = 1024;
p = 1024;

% dense matrix following normal distribution
A = randn(m,n);
B = randn(n,p);

% do projection based on the paper kyrillidis2014approximate
[PA, PB] = project(A,B, 'kyrillidis2014approximate',[1e-3, 1e-3, 1e-6]);

C_approx = PA*PB;


% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
C_norm = norm(C, 'fro');
AB_norm = A_norm * B_norm;
error_norm = norm(error, 'fro');

