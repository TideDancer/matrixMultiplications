% input A m*n, B n*p
k = log2(n);

% ------------- multiplication routing -------------




% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_nomr = norm(B, 'fro');
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
