% in Drines paper
% A m*n, B n*p
% sample c columns of A and corresponding rows of B

m = 1024;
n = 1024;
p = 1024;

A = randn(m,n);
B = randn(n,p);

% ------------- multiplication routing -------------
% the parameter l is to be tuned, here as default 1, but need to check paper
l = 1;
S = sampleElementL2(A, l);
R = sampleElementL2(B, l);

% sparsify representation significantly reduce computing time
S = sparse(S);
R = sparse(R);

C_approx = S*R;

% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
