% input A m*n, B n*p
n = 1024;

% dense matrix following normal distribution
A = randn(n,n);
B = randn(n,n);


% ------------- multiplication routing -------------
% build vector v
vec_v = randn(1,n);
V = zeros(n, n^2);
for i = 1:n
  V(i, n*i-n+1:n*i) = vec_v;
end
c = zeros(n^2, 1);

u = A*(B*v);

% solve the quation Vc = u
% using steepest descent method




% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_nomr = norm(B, 'fro');
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
